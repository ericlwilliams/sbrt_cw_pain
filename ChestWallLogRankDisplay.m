function ChestWallLogRankDisplay
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
%fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    
    
   
    
    %% Log-Rank p-value map
    
    
    % Vx
    f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
    
    pvx = CGobj{m}.mLogRank{f,2};% vmax
    f1 = pvx(:,:,6)== 0;
    f2 = pvx(:,:,5)<=0.05;
    pvx = pvx(:,:,5);
    pvx(~f1)=NaN;% Remove p-vals for which correlation not calc (unnecessary bc pval=1 otherwise)
    pvx(~f2)=NaN;
    
    fig_ctr=1;
    figure(fig_ctr);clf reset;
    fig_ctr=fig_ctr+1;
    % hack for better printing
    %gl = opengl('data')
    %opengl('software')
    
    h=imagesc(pvx');
    set(gca, 'ActivePositionProperty','OuterPosition');
    set(h,'alphadata',~isnan(pvx'));
    set(gca,'YDir','normal');
    
    
    xlim([0,60]);
    ylim([0,800]);
    colorbar;
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    set(gca,'xminortick','on','yminortick','on');
    
    xlabel('Dose [Gy]','fontsize',14);
    ylabel('Volume [cc]','fontsize',14);
    title('Significant Log-Rank p-values for Vx','fontsize',14);
    
    
    % VDx
    dose=39;
    
    %
    %     Vx=zeros(length(CGobj{m}.mBinsDose),CGobj{m}.mNumInGrp,1);
    %     for d=1:length(CGobj{m}.mBinsDose),
    %         for k=1:CGobj{m}.mNumInGrp
    %             Vx(d,k) = CGobj{m}.mGrp(k).fVolAtDose( d );
    %         end
    %     end
    %
    
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    % change to VDx after re-run, doesn't matter which one, all have same
    % info!
    f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};% vmax
    pvx = squeeze(pvx(fdose,:,:)); % look at volumes V_{fdose}, then pvx = [V_{fdose} bins]X[n1,c1,n2,c2,p,flg]
    f = pvx(:,6) == 0;
    
    % find minimum p-value and corresponding
    % volume split value  in 2-3rd quartile range
    v39_pvals = pvx(f,5);
    v39_vols = CGobj{m}.mBinsVol(f);
    
    [sorted_v39_pvals,idx_sorted_pvals] = sort(pvx(f,5));
    sorted_v39_vols = v39_vols(idx_sorted_pvals);
    
    mid_quart_pvals =...
        sorted_v39_pvals(find([sorted_v39_vols>quantile(sorted_v39_vols,0.25)].*[sorted_v39_vols<quantile(sorted_v39_vols,0.75)]));
    [min_mid_qrt_pval,~] = min(mid_quart_pvals);
    mid_qrt_best_split = find([v39_pvals==min_mid_qrt_pval]);
    if length(mid_qrt_best_split)>1
        mid_qrt_best_split = mid_qrt_best_split(1);
    end
    disp(['Best split v39 volume for 2-3rd quantiles: ',num2str(v39_vols(mid_qrt_best_split)),10,...
        ' with pval: ',num2str(min_mid_qrt_pval)]);
    
    
    
    %med_v39 = median(vx);
    figure(fig_ctr); clf reset;
    fig_ctr=fig_ctr+1;
    %[ax,~,~]=plotyy(CGobj{m}.mBinsVol(f),pvx(f,5),CGobj{m}.mBinsVol(f),pvx(f,7),@semilogy)
    semilogy(CGobj{m}.mBinsVol(f),pvx(f,5),'-','LineWidth',2); hold on;
    semilogy(CGobj{m}.mBinsVol(f),0.05,'r-','LineWidth',1);
    hold off; % grid on;
    xlabel('Volume (cc)','fontsize',14);
    ylabel('p-value','fontsize',14);
    title(strcat('Log-Rank p-values for V_{',num2str(dose),'} values'),'fontsize',12);
    [min_pvx,idx_pvx] = min(pvx(f,5));
    disp(['Minimum Log-Rank p-value is ',...
        num2str(min_pvx),...
        ', at V_{',num2str(dose),'} = ',...
        num2str(CGobj{m}.mBinsVol(idx_pvx))]);
    
    
    f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};% vmax
    min_pvx = nan(length(CGobj{m}.mBinsDose),1);
    idx_pvx = -ones(length(CGobj{m}.mBinsDose),1);
    med_vols = nan(length(CGobj{m}.mBinsDose),1);
    mn_vols = nan(length(CGobj{m}.mBinsDose),1);
    for i=1:length(CGobj{m}.mBinsDose)
        
        [~,fdose] = min(abs(CGobj{m}.mBinsDose - CGobj{m}.mBinsDose(i)));
        
        cur_pvx = squeeze(pvx(fdose,:,:)); % look at volumes V_{fdose}, then pvx = [V_{fdose} bins]X[n1,c1,n2,c2,p,flg]
        f = cur_pvx(:,6) == 0;
        if sum(f)<1
            continue;
        end
        
        Vx=zeros(CGobj{m}.mNumInGrp,1);
        for k=1:CGobj{m}.mNumInGrp
            Vx(k) = CGobj{m}.mGrp(k).fVolAtDose( CGobj{m}.mBinsDose(i));
        end
        
        
        %med_vols(i) = median(CGobj{m}.mBinsVol(f));
        %mn_vols(i) = mean(CGobj{m}.mBinsVol(f));
        med_vols(i) = median(Vx);
        mn_vols(i) = mean(Vx);
        
        cur_pvx =cur_pvx(f,5);
        [min_pvx(i),idx_pvx(i)] =min(cur_pvx);
    end
    
    figure(fig_ctr); clf reset;
    fig_ctr=fig_ctr+1;
    semilogy(CGobj{m}.mBinsDose,min_pvx,'-','LineWidth',2); hold on;
    semilogy(CGobj{m}.mBinsDose,0.05,'r-','LineWidth',1);
    set(gca,'xminortick','on','yminortick','on');
    grid on;
    xlabel('Dose [Gy]','fontsize',14);
    ylabel('Minimum Log-Rank p-value','fontsize',14);
    title('Log-Rank p-value vs V_{D}','fontsize',12);
    hold off; % grid on;
    
    f4=figure(fig_ctr); clf reset;
    fig_ctr=fig_ctr+1;
    %set(f4,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h=idx_pvx<1;
    idx_pvx(h)=[];
    hold on;
    j(1)=plot(CGobj{m}.mBinsDose(~h),CGobj{m}.mBinsVol(idx_pvx),'-','LineWidth',2);
    j(2)=plot(CGobj{m}.mBinsDose(~h),med_vols(~h),'r--','LineWidth',1);
    j(3)=plot(CGobj{m}.mBinsDose(~h),mn_vols(~h),'g--','LineWidth',1);
    hold off;
    
    set(gca,'xminortick','on','yminortick','on');
    grid on;
    legend(j,'Best V_{D}','Median V_{D}','Mean V_{D}');
    xlabel('Dose [Gy]','fontsize',14);
    ylabel('Volume of best Log-Rank split [cc]','fontsize',14);
    title('Log-Rank Volume split vs V_{D}','fontsize',12);
    hold off; % grid on;
    
    %% KM curve displays
    
    % compose complication KM curves
    
    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
    sa2=classKaplanMeierCurve(); % initialize a survivalanalysis obj
    CG = CGobj{m};
    % survival/complication time
    f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
    f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
    compdate = inf(CG.mNumInGrp,1);
    lastfollowup = inf(CG.mNumInGrp,1);
    compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
    lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
    compdate = min( lastfollowup, compdate );
    flgcensor = [CG.mGrp.mFlgCensor]';
    
    
    
    % V39
    %median split
    dose=39;
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};
    pvx = squeeze(pvx(fdose,:,:));
    f = pvx(:,6) < 2;
    
    
    vd=zeros(CG.mNumInGrp,1); % volume v at dose d
    numstart=CG.mLogRankMinSize;
    for d=fdose:fdose
        % volume under d
        vd(:)=0;
        for k=1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
        end
        g=find(vd>-1); % non-zeros volume cases
        
        % median volume
        vol = median(vd);
        %vol = 97;
        [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
        disp(['Median value of V_{',num2str(dose),'} = ',num2str(vol)]);
        disp(['Log-Rank p-value of V_{',num2str(dose),...
            '} split at ',num2str(vol),...
            ' = ',num2str(pvx(fvol,5))]);
        str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(fvol,5),'%10.2e\n'))};
        
        % (di,vj)
        numend=length(g)-numstart;
        for v=fvol:fvol
            % check smaple size
            flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
            if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                continue;
            end
            
            % assign properties of object sa
            survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
            fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
            sa.mSurvivalTime=survivedate;
            sa.mFlgCensor=fcensor;
            % compute survival curves and compare them
            sa=sa.fCalculateSurvivalCurve();
            sa=sa.fCombineSurvivalTime();
            sa=sa.fCompareSurvivalByLogrank();
        end
    end
    
    % plot KM curves
    figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    text(38,0.25,str_pval1,'FontSize',12);
    lgnd=legend(h_km,...
        strcat('V$_{39}\leq',num2str(vol,4),'$'),...
        strcat('V$_{39}\geq',num2str(vol,4),'$'),'Location','Best');
    set(lgnd,'FontSize',14);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('V_{39}, Median split','fontsize',14);
    
    dose=33;
    split=-1;
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    
    dose=39;
    split=-1;
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    
    dose=39;
    split=97;%best
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
%     % V33
%     %median split
%     dose=33;
%     [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
%     f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
%     pvx = CGobj{m}.mLogRank{f,2};
%     pvx = squeeze(pvx(fdose,:,:));
%     f = pvx(:,6) < 2;
%     
%     
%     vd=zeros(CG.mNumInGrp,1); % volume v at dose d
%     numstart=CG.mLogRankMinSize;
%     for d=fdose:fdose
%         % volume under d
%         vd(:)=0;
%         for k=1:CG.mNumInGrp
%             vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
%         end
%         g=find(vd>-1); % non-zeros volume cases
%         
%         % median volume
%         vol = median(vd);
%         [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
%         disp(['Median value of V_{',num2str(dose),'} = ',num2str(vol)]);
%         disp(['Log-Rank p-value of V_{',num2str(dose),...
%             '} split at ',num2str(vol),...
%             ' = ',num2str(pvx(fvol,5))]);
%         str_pval2 = {strcat('Log-Rank p-value = ',num2str(pvx(fvol,5),'%10.2e\n'))};
%         
%         % (di,vj)
%         numend=length(g)-numstart;
%         for v=fvol:fvol
%             % check smaple size
%             flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
%             if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
%                 continue;
%             end
%             
%             % assign properties of object sa
%             survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
%             fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();
%         end
%     end
%     
%     % plot KM curves
%     figure(6);  clf reset; hold on; % grid on;
%     h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
%     plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
%         1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
%     h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
%     plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
%         1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
%     %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
%     text(38,0.25,str_pval2,'FontSize',12);
%     lgnd=legend(h_km,...
%         strcat('V$_{33}\leq',num2str(vol,4),'$'),...
%         strcat('V$_{33}\geq',num2str(vol,4),'$'),'Location','Best');
%     set(lgnd,'FontSize',14);
%     h=legend;
%     set(h,'interpreter','latex');
%     
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel(['Months'],'fontsize',14);
%     ylabel(['Probability of CW Pain'],'fontsize',14);
%     title('V_{33}, Median split','fontsize',14);
%     
    
    % V30
    %median split
%     dose=30;
dose=30;
    
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
    cur_pvx = CGobj{m}.mLogRank{f,2};
    cur_pvx = squeeze(cur_pvx(fdose,:,:));
   f = cur_pvx(:,6) < 2;
    
      
    f7=figure(fig_ctr); clf reset;
    fig_ctr=fig_ctr+1;

    set(f7,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    x_axis = CGobj{m}.mBinsVol(f);
    x_axis = x_axis(1:150);
    v30_hrs = cur_pvx(f,7);
    v30_hrs = v30_hrs(1:150);
    v30_pvals = cur_pvx(f,5);
    v30_pvals = v30_pvals(1:150);
    %semilogy(CGobj{m}.mBinsVol(f),cur_pvx(f,5),'-','LineWidth',2); hold on;
    [ax,h1,h2]=plotyy(x_axis,v30_hrs,x_axis,v30_pvals,@semilogy);hold on;
    set(h1,'LineWidth',2);
    set(h2,'LineWidth',2);
    set(ax(1),'YLim',[2.5 7]);
    set(ax(1),'YTick',[2.5:1:7]);
    set(get(ax(1),'Ylabel'),'String','Hazard Ratio','FontSize',15);
        
    set(get(ax(2),'Ylabel'),'String','p-value','FontSize',15);
    semilogy(x_axis,0.05,'r-','LineWidth',1);
    %xlim([0 max(x_axis)]);
    set(ax(1),'XLim',[0 max(x_axis)]);
    set(ax(2),'XLim',[0 max(x_axis)]);
    set(ax,'FontSize',14);
    %line(vol,ylim);
    %here
    Lmed = line([med_vols(31) med_vols(31)],ylim,'Color','k','LineWidth',2);
    %L30 = line([30 30],ylim,'Color','g','LineWidth',2);
    %L50 = line([50 50],ylim,'Color','c','LineWidth',2);
    %L70 = line([70 70],ylim,'Color','m','LineWidth',2);
    hold off; % grid on;
    xlabel('Volume (cc)','fontsize',15);
    %ylabel('p-value','fontsize',14);
    title(strcat('Log-Rank p-values for V_{',num2str(dose),'} values'),'fontsize',14);
%    lgnd=legend([Lmed L30 L50 L70],['Median = ',num2str(med_vols(30))],'30cc','50cc','70cc','Location','SouthEast');
    lgnd=legend(Lmed,['Median = ',num2str(med_vols(31),4),' cc']);
    set(lgnd,'FontSize',14);
    
%     vd=zeros(CG.mNumInGrp,1); % volume v at dose d
%     numstart=CG.mLogRankMinSize;
%     for d=fdose:fdose
%         % volume under d
%         vd(:)=0;
%         for k=1:CG.mNumInGrp
%             vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
%         end
%         g=find(vd>-1); % non-zeros volume cases
%         
%         % median volume
%         vol = median(vd);
%         %vol = 50;
%         [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
%         disp(['Median value of V_{',num2str(dose),'} = ',num2str(vol)]);
%         disp(['Log-Rank p-value of V_{',num2str(dose),...
%             '} split at ',num2str(vol),...
%             ' = ',num2str(cur_pvx(fvol,5))]);
%         str_pval2 = {strcat('Log-Rank p-value = ',num2str(cur_pvx(fvol,5),'%10.2e\n'))};
%         
%         
%         
%         % (di,vj)
%         numend=length(g)-numstart;
%         for v=fvol:fvol
%             % check smaple size
%             flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
%             if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
%                 continue;
%             end
%             
%             % assign properties of object sa
%             survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
%             fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();
%         end
%         
%         
%         vol2 = 30; %cc
%         [~,fvol2] = min(abs(CGobj{m}.mBinsVol - vol2));
%         disp(['Median value of V_{',num2str(dose),'} = ',num2str(vol2)]);
%         disp(['Log-Rank p-value of V_{',num2str(dose),...
%             '} split at ',num2str(vol2),...
%             ' = ',num2str(cur_pvx(fvol2,5))]);
%         str_pval3 = {strcat('Log-Rank p-value = ',num2str(cur_pvx(fvol2,5),'%10.2e\n'))};
%         
%         for v=fvol2:fvol2
%             % check smaple size
%             flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
%             if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
%                 continue;
%             end
%             
%             % assign properties of object sa
%             survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
%             fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
%             sa2.mSurvivalTime=survivedate;
%             sa2.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa2=sa2.fCalculateSurvivalCurve();
%             sa2=sa2.fCombineSurvivalTime();
%             sa2=sa2.fCompareSurvivalByLogrank();
%         end
%         
%         
%         
%     end
%     
%     
%     % plot KM curves
%     % V_{30} split at median: 50cc
%     f8=figure(8);  clf reset; hold on; % grid on;
%     set(f8,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%     h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
%     plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
%         1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
%     h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
%     plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
%         1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
%     %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
%     text(38,0.25,str_pval2,'FontSize',14);
%     ylim([0 0.53]);
%     lgnd=legend(h_km,...
%         strcat('V$_{30}\leq',num2str(vol,4),'$~cc'),...
%         strcat('V$_{30}\geq',num2str(vol,4),'$~cc'),'Location','Best');
%     set(lgnd,'FontSize',14);
%     h=legend;
%     set(h,'interpreter','latex');
%     
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel(['Months'],'fontsize',14);
%     ylabel(['Probability of CW Pain'],'fontsize',14);
%     title(strcat('V_{30}, Median split (',num2str(vol),' cc)'),'fontsize',14);
%     
%     % plot KM curves
%     % V_{30} split at 30cc
%     f9=figure(9);  clf reset; hold on; % grid on;
%     set(f9,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%     
%     h_km(1)=stairs(sa2.mSurvivalTimeSorted{1},1-sa2.mSurvivalCurve{1});
%     plot(sa2.mSurvivalTimeSorted{1}(sa2.mCensorStatistics{1}(:,1)),...
%         1-sa2.mSurvivalCurve{1}(sa2.mCensorStatistics{1}(:,1)),'+');
%     h_km(2)=stairs(sa2.mSurvivalTimeSorted{2},1-sa2.mSurvivalCurve{2},'r');
%     plot(sa2.mSurvivalTimeSorted{2}(sa2.mCensorStatistics{2}(:,1)),...
%         1-sa2.mSurvivalCurve{2}(sa2.mCensorStatistics{2}(:,1)),'r+');
%     %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
%     text(38,0.25,str_pval3,'FontSize',14);
%     lgnd=legend(h_km,...
%         strcat('V$_{30}\leq',num2str(vol2,4),'$~cc'),...
%         strcat('V$_{30}\geq',num2str(vol2,4),'$~cc'),'Location','Best');
%     set(lgnd,'FontSize',14);
%     h=legend;
%     set(h,'interpreter','latex');
%     ylim([0 0.53]);
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel(['Months'],'fontsize',14);
%     ylabel(['Probability of CW Pain'],'fontsize',14);
%     title('V_{30}, 30cc split','fontsize',14);
    
    dose=30;
    split=med_vols(31);
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    ylim([0 0.65]);
    
    dose=30;
    split=50;
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    ylim([0 0.65]);
    
    
    dose=30;
    split=30;
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    ylim([0 0.53]);
    
    dose=30;
    split=70;
    fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
    fig_ctr=fig_ctr+1;
    ylim([0 0.53]);
    
    
    %% BMI split
% at median
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 
    bmi_data = bmi_data(bmi_idx);
    
    numstart=CG.mLogRankMinSize;
            %here
    
    bmi_split = median(bmi_data);
   
     
     
        % (di,vj)
    numend=length(bmi_data)-numstart;
        
    flg_volbelow1= (bmi_data<bmi_split);
    f=length(find(flg_volbelow1)); % group DVHs by (d,v)
    if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
        continue;
    end
    
    % assign properties of object sa
%     survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
%     fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
%     sa.mSurvivalTime=survivedate;
%     sa.mFlgCensor=fcensor;
%     % compute survival curves and compare them
%     sa=sa.fCalculateSurvivalCurve();
%     sa=sa.fCombineSurvivalTime();
%     sa=sa.fCompareSurvivalByLogrank();
    
    [~,fsplit] = min(abs(unique(bmi_data)-bmi_split));
    f = cellfun(@(x) strcmpi('BMI',x),CGobj{m}.mLogRank(:,1));
    if ~isempty(f)
    pvx = CGobj{m}.mLogRank{f,2};
    
    %here
    sa = pvx{2}{fsplit};
    
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    
    
    % plot KM curves
    f11=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f11,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    ylim([0 0.55]);
    text(38,0.12,str_pval1,'FontSize',16);
    lgnd=legend(h_km,...
        strcat('$\rm{BMI} <',num2str(bmi_split,4),'$'),...
        strcat('$\rm{BMI} \geq',num2str(bmi_split,4),'$'),...
        'Location','Best');
    set(lgnd,'FontSize',16);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('BMI, Median split','fontsize',14);
    
    f = cellfun(@(x) strcmpi('BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};% vmax
    pvx = pvx{1}
    %pvx = squeeze(pvx(fdose,:,:)); % look at volumes V_{fdose}, then pvx = [V_{fdose} bins]X[n1,c1,n2,c2,p,flg]
    f = pvx(:,6) == 0;
    f12=figure(fig_ctr); clf reset;
    fig_ctr=fig_ctr+1;
    set(f12,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    unique_bmi_data = unique(bmi_data);
    [min_pvx,idx_pvx] = min(pvx(f,5));
    unique_bmi_data = unique_bmi_data(f);
    semilogy(unique_bmi_data,pvx(f,5),'-','LineWidth',2); hold on;
    semilogy([min(unique_bmi_data) max(unique_bmi_data)],[0.05 0.05],'r--','LineWidth',2);
    semilogy([unique_bmi_data(idx_pvx) unique_bmi_data(idx_pvx)],ylim,'g--','LineWidth',2);
    hold off; % grid on;
    xlim([min(unique_bmi_data) max(unique_bmi_data)]);
    xlabel('BMI','fontsize',14);
    ylabel('p-value','fontsize',14);
    title(strcat('Log-Rank p-values for BMI split values'),'fontsize',14);
    
    str_pval1 = {strcat('Log-Rank p-value = ',...
                        num2str(min_pvx,'%10.2e\n'),...
                        ' at BMI = ',...
                        num2str(unique_bmi_data(idx_pvx)))};
    text(15,0.001,str_pval1,'FontSize',14);      

    %% BMI at min pval
    best_bmi = unique_bmi_data(idx_pvx);
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 
    bmi_data = bmi_data(bmi_idx);
    
    numstart=CG.mLogRankMinSize;
            %here
    
    bmi_split = best_bmi;
        
        % (di,vj)
    numend=length(bmi_data)-numstart;
        
    flg_volbelow1= (bmi_data<bmi_split);
    f=length(find(flg_volbelow1)); % group DVHs by (d,v)
    if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
        continue;
    end
    
    cox_beta=coxphfit(~flg_volbelow1,compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    cox_hr = exp(cox_beta);
    disp(['HR(BMI,min p-val): ',num2str(cox_hr)]);
    
    % assign properties of object sa
%     survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
%     fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
%     sa.mSurvivalTime=survivedate;
%     sa.mFlgCensor=fcensor;
%     % compute survival curves and compare them
%     sa=sa.fCalculateSurvivalCurve();
%     sa=sa.fCombineSurvivalTime();
%     sa=sa.fCompareSurvivalByLogrank();
%     
    [~,fsplit] = min(abs(unique(bmi_data)-bmi_split));
    f = cellfun(@(x) strcmpi('BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};
    % here
    sa = pvx{2}{fsplit};
    
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    
    
    % plot KM curves
    f13=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f13,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    ylim([0 0.55]);
    text(38,0.15,str_pval1,'FontSize',16);
    lgnd=legend(h_km,...
        strcat('$\rm{BMI} <',num2str(bmi_split,4),'$'),...
        strcat('$\rm{BMI} \geq',num2str(bmi_split,4),'$'),...
        'Location','Best');
    set(lgnd,'FontSize',16);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('BMI, best split','fontsize',14);
    
      %% BMI at 29
    best_bmi = unique_bmi_data(idx_pvx);
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 
    bmi_data = bmi_data(bmi_idx);
    
    numstart=CG.mLogRankMinSize;
            %here
    
    %bmi_split = 29;
    bmi_split = 30;
        
        % (di,vj)
    numend=length(bmi_data)-numstart;
        
    flg_volbelow1= (bmi_data<bmi_split);
    f=length(find(flg_volbelow1)); % group DVHs by (d,v)
    if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
        continue;
    end
    
    % assign properties of object sa
%     survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
%     fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
%     sa.mSurvivalTime=survivedate;
%     sa.mFlgCensor=fcensor;
%     % compute survival curves and compare them
%     sa=sa.fCalculateSurvivalCurve();
%     sa=sa.fCombineSurvivalTime();
%     sa=sa.fCompareSurvivalByLogrank();
    
    [~,fsplit] = min(abs(unique(bmi_data)-bmi_split));
    f = cellfun(@(x) strcmpi('BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};
    
    %here
    sa = pvx{2}{fsplit};
    cox_beta=coxphfit(~flg_volbelow1,compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    cox_hr = exp(cox_beta);
    disp(['Cox HR for BMI: ',num2str(cox_hr)]);
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    %str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    str_pval1 = ['Log-Rank p-value = ',num2str(pvx(1,5),'%3.1e\n'),10,...
    'HR = ',num2str(cox_hr,'%3.1f')];
    % plot KM curves
    f14=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f14,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    ylim([0 0.65]);
    textbp(str_pval1,'FontSize',16);
    lgnd=legend(h_km,...
        strcat('$\rm{BMI} <',num2str(bmi_split,4),'$'),...
        strcat('$\rm{BMI} \geq',num2str(bmi_split,4),'$'),...
        'Location','Best');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'fontsize',16);
    xlabel(['Years'],'fontsize',20);
    ylabel(['Incidence of CW Pain'],'fontsize',20);
    title(['BMI, split at ',num2str(bmi_split,4)],'fontsize',20);
    
    %% V_{39} + BMI split


    dose=39;         
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 

    vd=zeros(CG.mNumInGrp,1); % volume v at dose d
    numstart=CG.mLogRankMinSize;
    for d=fdose:fdose
        % volume under d
        vd(:)=0;
        for k=1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
        end
        %g=find(vd>-1); % non-zeros volume cases
        g=find(bmi_idx);
        
        %here
        vd_bmi = 0.0207.*vd(g) + 0.0408.*bmi_data(g);
        vd_bmi_split = median(vd_bmi);
        
        
        %disp(['Median value of V_{',num2str(dose),'} = ',num2str(vol)]);
        %disp(['Log-Rank p-value of V_{',num2str(dose),...
         %   '} split at ',num2str(vol),...
          %  ' = ',num2str(pvx(fvol,5))]);

        
        
        % (di,vj)
        numend=length(g)-numstart;
        
        flg_volbelow1= (vd_bmi<vd_bmi_split);
        f=length(find(flg_volbelow1)); % group DVHs by (d,v)
        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
            continue;
        end
%             
%             % assign properties of object sa
%             survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
%             fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();

    end
            
    [~,fsplit] = min(abs(unique(vd_bmi)-vd_bmi_split));
    f = cellfun(@(x) strcmpi('V39BMI',x),CGobj{m}.mLogRank(:,1));
    
    if ~isempty(f)
    pvx = CGobj{m}.mLogRank{f,2};

    sa = pvx{2}{fsplit};
    
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    end
    
    % plot KM curves
    f10=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f10,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    text(38,0.25,str_pval1,'FontSize',16);
    ylim([0 0.8]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{39}}\times V_{39}+\beta_{BMI}\times BMI <',num2str(vd_bmi_split,3),'$'),...
        strcat('$\beta_{V_{39}}\times V_{39}+\beta_{BMI}\times BMI \geq',num2str(vd_bmi_split,3),'$'),...
        'Location','NorthEast');
    set(lgnd,'FontSize',14);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('V_{39} + BMI, Median split','fontsize',14);
    end    
        %% V_{30} + BMI split


    dose=30;         
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 

    vd=zeros(CG.mNumInGrp,1); % volume v at dose d
    numstart=CG.mLogRankMinSize;
    for d=fdose:fdose
        % volume under d
        vd(:)=0;
        for k=1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
        end
        %g=find(vd>-1); % non-zeros volume cases
        g=find(bmi_idx);
        
        %here
        vd_bmi = 0.0131*vd(g) + 0.0426*bmi_data(g);
        vd_bmi_split = median(vd_bmi);
        
        numend=length(g)-numstart;
        
        flg_volbelow1= (vd_bmi<vd_bmi_split);
        f=length(find(flg_volbelow1)); % group DVHs by (d,v)
        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
            continue;
        end

%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();
            
            
    [~,fsplit] = min(abs(unique(vd_bmi)-vd_bmi_split));
    f = cellfun(@(x) strcmpi('V30BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};

    sa = pvx{2}{fsplit};
    
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    end
    
    % plot KM curves
    f16=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f16,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    %text(38,0.25,str_pval1,'FontSize',16);
    textbp(str_pval1,'FontSize',16);
    ylim([0 0.65]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI <',num2str(vd_bmi_split,3),'$'),...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \geq',num2str(vd_bmi_split,3),'$'),...
        'Location','NorthEast');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'fontsize',16);
    xlabel(['Months'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    title('V_{30} + BMI, Median split','fontsize',20);
    
    %% V_{30} + BMI + Tx

    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj

    dose=30;
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0; 

     [txCox,~,~] = CGobj{m}.fCoxParameter_DVH('Tx');
    tx_data = [txCox.data_exposure];

    
    vd=zeros(CG.mNumInGrp,1); % volume v at dose d
    numstart=CG.mLogRankMinSize;
    for d=fdose:fdose
        % volume under d
        vd(:)=0;
        for k=1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
        end
        %g=find(vd>-1); % non-zeros volume cases
        g=find(bmi_idx);
        
        %here
        vd_bmi_tx = 0.011902*vd(g) +...
            0.042235*bmi_data(g) +...
            0.000527*tx_data(g);
        
        vd_bmi_tx_split = median(vd_bmi_tx);

            
             
        
        % (di,vj)
        numend=length(g)-numstart;
        
        flg_volbelow1= (vd_bmi_tx<vd_bmi_tx_split);
        f=length(find(flg_volbelow1)); % group DVHs by (d,v)
        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
            continue;
        end
            
            % assign properties of object sa
%             survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
%             fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();
            
            
    [~,fsplit] = min(abs(unique(vd_bmi_tx)-vd_bmi_tx_split));
    f = cellfun(@(x) strcmpi('V30BMITX',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};

    
    sa = pvx{2}{fsplit};

    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(pvx(1,5),'%10.2e\n'))};
    
    end
    
    % plot KM curves
    f15=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f15,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    %text(38,0.25,str_pval1,'FontSize',16);
    
    ylim([0 0.65]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI+\beta_{Tx}\times Tx<',num2str(vd_bmi_tx_split,3),'$'),...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI+\beta_{Tx}\times Tx \geq',num2str(vd_bmi_tx_split,3),'$'),...
        'Location','NorthEast');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    textbp(str_pval1,'FontSize',18);
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'fontsize',16);
    xlabel(['Years'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    title('V_{30} + BMI + Tx, Median split','fontsize',20);
    
    
    %% V30 + BMI + TX best split
   f = cellfun(@(x) strcmpi('V30BMITX',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};

    [min_pval,min_pvx_idx] = min(pvx{1}(:,5));
    
    sa = pvx{2}{min_pvx_idx};
    vd_bmi_tx = 0.011902*vd(bmi_idx) +...
                0.042235*bmi_data(bmi_idx) +...
                0.000527*tx_data(bmi_idx);
    vd_bmi_tx_split = vd_bmi_tx(min_pvx_idx);
            
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(min_pval,'%10.2e\n'))};
    

    
    % plot KM curves
    f18=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f18,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    textbp(str_pval1,'FontSize',16);
    ylim([0 0.8]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI+\beta_{Tx}\times Tx<',num2str(vd_bmi_tx_split,3),'$'),...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI+\beta_{Tx}\times Tx \geq',num2str(vd_bmi_tx_split,3),'$'),...
        'Location','East');
    set(lgnd,'FontSize',14);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('V_{30} + BMI + Tx, best split','fontsize',14);
    
     %% V30 + BMI  best split
   f = cellfun(@(x) strcmpi('V30BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};

    [min_pval,min_pvx_idx] = min(pvx{1}(:,5));
    
    sa = pvx{2}{min_pvx_idx};
    vd_bmi = 0.0131*vd(bmi_idx) + 0.0426*bmi_data(bmi_idx);
    vd_bmi_split = vd_bmi(min_pvx_idx);
            
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(min_pval,'%10.2e\n'))};
    

    
    % plot KM curves
    f19=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f19,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    ylim([0 0.65]);
    textbp(str_pval1,'FontSize',16);
    %ylim([0 0.8]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI <',num2str(vd_bmi_split,3),'$'),...
        strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \geq',num2str(vd_bmi_split,3),'$'),...
        'Location','Best');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    title('V_{30} + BMI, best split','fontsize',20);
    
      %% V39 + BMI  best split
     f = cellfun(@(x) strcmpi('V39BMI',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};

    [min_pval,min_pvx_idx] = min(pvx{1}(:,5));
    
    sa = pvx{2}{min_pvx_idx};
    vd_bmi = 0.0207.*vd(bmi_idx) + 0.0408.*bmi_data(bmi_idx);
    vd_bmi_split = vd_bmi(min_pvx_idx);
            
    pvx = pvx{1};
    pvx = squeeze(pvx(fsplit,:,:));
    f = pvx(:,6) < 2;
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(min_pval,'%10.2e\n'))};
    

    
    % plot KM curves
    f20=figure(fig_ctr);  clf reset; hold on; % grid on;
    fig_ctr=fig_ctr+1;
    set(f20,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    text(38,0.1,str_pval1,'FontSize',16);
    ylim([0 0.8]);
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{39}}\times V_{39}+\beta_{BMI}\times BMI <',num2str(vd_bmi_split,3),'$'),...
        strcat('$\beta_{V_{39}}\times V_{39}+\beta_{BMI}\times BMI \geq',num2str(vd_bmi_split,3),'$'),...
        'Location','East');
    set(lgnd,'FontSize',14);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title('V_{39} + BMI, best split','fontsize',14);
end

    
toc;
end