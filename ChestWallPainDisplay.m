function ChestWallPainDisplay
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\Desktop\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';

if isunix
    fp=strrep(fp,'G:','/media/SKI_G');
end
fn = {%'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
      %'RIMNER_BC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
      %'RIMNER_AC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
      'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
      %'MUTTER_BC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
      %'MUTTER_AC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};

      %test
      fn={'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
          
CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end


for m = 1:length(fn)
    disp(fn{m});
    CGobj{m}.fCoxFig_DVH('VDx');
    set(figure(1),'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(figure(1),'Name',fn{m});
    set(figure(2),'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(figure(2),'Name',fn{m});
    %CGobj{m}.fCoxFig_DVH('Vx');
    continue;
    
    figure(4); clf reset;
    CGobj{m}.fCoxRiskVDxFig_DVH();
    xlabel(''); ylabel('');
    figure(5);
    CGobj{m}.fCoxRiskVDxFig_DVH(30);
    xlabel(''); ylabel('');
    
end




% plot DVH groups
disp('DVH groups');
for m = 1:length(fn)
    cur_f=figure(m); clf reset; hold on;
    set(cur_f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(cur_f,'Name',fn{m});
    CGobj{m}.fDVHCurvesSummary_DVH();
    xlabel('Dose [Gy]','FontSize',18);ylabel('Vol [cc]','FontSize',18);
end





% plot DVHs
disp('DVH curves:');
for m = 1:length(fn)
    disp(fn{m});
    f1=figure(m); clf reset; hold on; % grid on;
    set(f1,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    f = [CGobj{m}.mGrp.mFlgCensor];
    % DVHs of censored patients
    g = find(f);
    for k = 1:length(g)
        plot(CGobj{m}.mGrp(g(k)).mDoseBins_org, CGobj{m}.mGrp(g(k)).mVolCum);
    end
    % DVHs of complicated patients
    g = find(~f);
    for k = 1:length(g)
        plot(CGobj{m}.mGrp(g(k)).mDoseBins_org, CGobj{m}.mGrp(g(k)).mVolCum,'r');
    end
    xlabel('Dose [Gy]','FontSize',18);
    ylabel('Vol [cc]','FontSize',18);
    set(gca,'xminortick','on','yminortick','on');

end




% median complication time
%mediantime(cw_dataset,cur_comp);



% plot V22Gy vs. complication
disp('CW pain vs. V22Gy');
m = 1;
d = 22;
% extract V22Gy
vd=zeros(CGobj{m}.mNumInGrp,1); % volume v at dose d
% volume under d
for k=1:CGobj{m}.mNumInGrp
    vd(k) = CGobj{m}.mGrp(k).fVolAtDose(d);
end
% complication
flgcensor = [CGobj{m}.mGrp.mFlgCensor]';
% plot
figure(m); clf reset;
plot(vd(flgcensor),0,'b*', vd(~flgcensor),1,'r*');

% survival curves
disp('Survival curves');
m=1;
figure(m+m-1); clf reset; hold on; % grid on;
sa = CGobj{m}.mKaplanMeierSurvivalOverall;
stairs(sa.mSurvivalTimeSorted{1},sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
set(gca,'Ylim',[0,1]);
%         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
set(gca,'xminortick','on','yminortick','on');
xlabel(['Months']);
ylabel(['Overall survival probability']);
% median survival time
disp(['median survival time: ',num2str(median(sa.mSurvivalTime{1}))]);


% complication incidence curve
disp('Complication incidence curves');
m = 1;
cur_fig=figure(m+100); clf reset; hold on; % grid on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

sa = CGobj{m}.mKaplanMeierCompOverall;
stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'LineWidth',2);
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',15);
%         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
set(gca,'xminortick','on','yminortick','on');
xlabel('Years','FontSize',18);
ylabel('Inicidence of grade >= 2 Chestwall Pain','FontSize',18);
set(gca,'FontSize',16);
text(3,0.1,['Median onset time: ',num2str(median(sa.mSurvivalTime{1}(~sa.mFlgCensor{1}))/12,2),' yr'],'FontSize',18,'BackgroundColor','w');
grid on;

% median incident time
disp(['meidan incident time: ', num2str(median(sa.mSurvivalTime{1}(~sa.mFlgCensor{1})))]);


set(cur_fig,'Color','w');
export_fig(cur_fig,...
    [fig_loc,'cwp_cuminc'],'-pdf');
disp(['Saving ',fig_loc,'cwp_cuminc.pdf']);


        
% Relapse free
% figure(m+m); clf reset; hold on; % grid on;
% sa = CGobj{m}.mKaplanMeierRelapseFree;
% stairs(sa.mSurvivalTimeSorted{1},sa.mSurvivalCurve{1});
% plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
%     sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
% set(gca,'Ylim',[0,1]);
% %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
% set(gca,'xminortick','on','yminortick','on');


% Cox model p-value curves
disp('Cox model p-values:');
p_significant=0.05;
for m = 1:length(fn)
    disp(fn{m});
    % VDx
    disp('VDx');
    figure(1); clf reset;
    [pCox,fCox] = CoxP(CGobj{m},'VDx');
    semilogy(CGobj{m}.mBinsDose,pCox,'b.-'); hold on;
    %         semilogy(CGobj{m}.mBinsDose(fCox),pCox(fCox),'r*');
    semilogy(CGobj{m}.mBinsDose,repmat(p_significant,size(CGobj{m}.mBinsDose)),'r--');
    hold off; % grid on;
    set(gca,'xminortick','on','yminortick','on');
    disp(['p-values: V20: ',num2str(pCox(21)),'  V30: ', num2str(pCox(31))]);
    % DVx
    disp('DVx');
    figure(2); clf reset;
    [pCox,fCox] = CoxP(CGobj{m},'DVx');
    semilogy(CGobj{m}.mBinsVol,pCox,'b.-'); hold on;
    semilogy(CGobj{m}.mBinsVol,repmat(p_significant,size(CGobj{m}.mBinsVol)),'r--');
    hold off; % grid on;
    set(gca,'xminortick','on','yminortick','on');
    
    
    
    % confidence interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the following code heavily depends on the curve's shape and has to be modified accordingly%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('CI');
    figure(3); clf reset;
    [allCox,flgCox,~] = CGobj{m}.fCoxParameter_DVH('VDx');
    logl = [allCox.logl];
    f = [allCox(flgCox).p]>0.3;
    g = find(flgCox);
    logl(g(f)) = -inf; % modify the loglikelihood of apparently non-significant points to avoid their disturbance in the following computation
    [mx,mxloc] = max(logl); % the maximum log likelihood and its location
    pct68 = mx - 0.5; % 68% confidence
    pct95 = mx - 1.96; % 95% confidence
    plot(CGobj{m}.mBinsDose(flgCox),[allCox(flgCox).logl],'b.-'); hold on;
    plot(CGobj{m}.mBinsDose(flgCox),repmat(pct68,size(CGobj{m}.mBinsDose(flgCox))),'r--');
    plot(CGobj{m}.mBinsDose(flgCox),repmat(pct95,size(CGobj{m}.mBinsDose(flgCox))),'c--');
    hold off;
    set(gca,'xminortick','on','yminortick','on');
    % compute the CI intervals
    % 68%
    % left end point
    g = find([logl(1:mxloc)] <= pct68);
    if isempty(g) % the curve is above pct68
        pct68l = -1;
    else % the curve climbed from below to above pct68
        fg = [g(end) g(end)+1];
        pct68l = interp1([logl(fg)],CGobj{m}.mBinsDose(fg),pct68);
    end
    % right end point
    g = [logl(mxloc+1:end)] >= pct68;
    if all(g) % the curve is above pct68
        pct68r = -1;
    else % the curve slide down below pct68
        g = find(g)+mxloc;
        fg = [g(end) g(end)+1];
        pct68r = interp1([logl(fg)],CGobj{m}.mBinsDose(fg),pct68);
    end
    disp(['68% interval: ', num2str([pct68l, pct68r])]);
    
    % 95%
    % left end point
    g = find([logl(1:mxloc)] <= pct95);
    if isempty(g) % the curve is above pct95
        pct95l = -1;
    else % the curve climbed from below to above pct95
        fg = [g(end) g(end)+1];
        pct95l = interp1([logl(fg)],CGobj{m}.mBinsDose(fg),pct95);
    end
    % right end point
    g = [logl(mxloc+1:end)] >= pct95;
    if all(g) % the curve is above pct95
        pct95r = -1;
    else % the curve slide down below pct95
        g = find(g)+mxloc;
        fg = [g(end) g(end)+1];
        pct95r = interp1([logl(fg)],CGobj{m}.mBinsDose(fg),pct95);
    end
    disp(['95% interval: ', num2str([pct95l, pct95r])]);
    
    % h(t) for V20 and V30
    % V20
    figure(4); clf reset; hold on;
    d = 20;
    [~,f] = min(abs(CGobj{m}.mBinsDose - d));
    plot(allCox(f).h(:,1),allCox(f).h(:,2),'b');
    % V30
    d = 30;
    [~,f] = min(abs(CGobj{m}.mBinsDose - d));
    plot(allCox(f).h(:,1),allCox(f).h(:,2),'r');
    set(gca,'xminortick','on','yminortick','on');
    
    % response function for V20
end

% % check p-values at 30cc
%     dose = 30;
%     vol = 30;
%     for m = 1:length(fn)

% KM curve for (30Gy 71cc)
% search the p-values for Vx with x = 30 Gy
%dose = 30;
%vol = 71;
dose = 40;
vol = 30;
for m = 1:length(fn)
    [~,fdose] = min(abs(CGobj{m}.mBinsDose - dose));
    f = cellfun(@(x) strcmpi('DVx',x),CGobj{m}.mLogRank(:,1));
    pvx = CGobj{m}.mLogRank{f,2};
    pvx = squeeze(pvx(fdose,:,:));
    f = pvx(:,6) < 2;
    figure(m); clf reset;
    semilogy(CGobj{m}.mBinsVol(f),pvx(f,5),'-','LineWidth',2); hold on;
    semilogy(CGobj{m}.mBinsVol(f),0.05,'r-','LineWidth',1);
    hold off; % grid on;
    set(gca,'xminortick','on','yminortick','on');
    xlabel('Volume (cc)'); ylabel('p-value');
    %     [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
    %     disp(CGobj{m}.mBinsVol(fvol));
    %     disp(pvx(fvol:fvol+5,:));
    
    % compose complication KM curves
    [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
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
    % Vx
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
        %                     vol = median(vd);
        [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
        disp(CGobj{m}.mBinsVol(fvol));
        disp(pvx(fvol-2:fvol+2,:));
        
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
    figure(m);  clf reset; hold on; % grid on;
    stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months']);
    ylabel(['Probability of CW Pain']);
end

% p-values of medium line log-rank test
for m = 1:length(fn)
    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
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
    
    dosebins = 0:max(CG.mBinsDose);
    pvalues = zeros(length(dosebins),1);
    vd = zeros(CG.mNumInGrp,1);
    for n = 1:length(dosebins)
        vd(:) = 0;
        for k = 1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose(dosebins(n));
        end
        vm = median(vd);
        flg_smallvol = vd<vm;
        survivedate = {compdate(flg_smallvol); compdate(~flg_smallvol)};
        fcensor = {flgcensor(flg_smallvol); flgcensor(~flg_smallvol)};
        sa.mSurvivalTime=survivedate;
        sa.mFlgCensor=fcensor;
        try
            sa=sa.fCalculateSurvivalCurve();
            sa=sa.fCombineSurvivalTime();
            sa=sa.fCompareSurvivalByLogrank();
            pvalues(n) = sa.mpValue;
        catch
        end
    end
    % plot the p-values
    figure(m); clf reset;
    semilogy(dosebins,pvalues,'b-'); hold on;
    semilogy(dosebins,0.05,'r-','LineWidth',1); hold off;
    set(gca,'xminortick','on','yminortick','on');
end
toc;
end


function [pCox,flganti] = CoxP(CGobj,strCoxVx)
f = cellfun(@(x) strcmpi(strCoxVx,x),CGobj.mCoxParameter(:,1)); % search the label
if any(f)
    % p-values in Cox model
    pCox = [CGobj.mCoxParameter{f,2}.logl]'; % some logl is infinite, indicating no data for those values
    fCox = isfinite(pCox);
    pCox(fCox) = [CGobj.mCoxParameter{f,2}(fCox).p]';
    pCox(~fCox)=1;
    % correlation
    f = [CGobj.mCoxParameter{f,2}.beta]';
    flganti = f<0;
else
    pCox = []; flganti = [];
end
end


function mediantime(fn,cur_comp)
% patient raw data
%fn='Z:\elw\MATLAB\original_data\CW_old\20100424CWFinalCohort';

if isunix
    fn=strrep(fn,'G:','/media/SKI_G');
end
load(fn,'xlsraw');
PtInfo = classDataFromXls();
PtInfo.mXlsRaw = xlsraw;

% baseline
PtInfo.mLabel = 'IGRT End Date';
PtInfo = PtInfo.fExtractColData();
flgdateIGRT = PtInfo.mFlg;
dateIGRT = zeros(size(flgdateIGRT));
f = PtInfo.mData(flgdateIGRT);
dateIGRT(flgdateIGRT) = datenum(f);

% grade 2
% complication grade 2
PtInfo.mLabel = cur_comp;
PtInfo = PtInfo.fExtractColData();
flgpain2 = PtInfo.mFlg;
datepain2 = inf(size(flgpain2));
datepain2(flgpain2) = datenum(PtInfo.mData(flgpain2));

% grade 3
PtInfo.mLabel = 'Pain Grade 3';
PtInfo = PtInfo.fExtractColData();
flgpain3 = PtInfo.mFlg;
datepain3 = inf(size(PtInfo.mData));
datepain3(flgpain3) = datenum(PtInfo.mData(flgpain3));

% first onset of pain
datepain = (min(datepain2,datepain3) - dateIGRT)/30;
f = isfinite(datepain);
disp(['median, min, and max of first onset pain: ', num2str([median(datepain(f)), min(datepain(f)), max(datepain(f))])]);

% maxium onset of pain
datepain = datepain2;
datepain(flgpain3) = datepain3(flgpain3);
datepain = (datepain - dateIGRT)/30;
f = isfinite(datepain);
disp(['median, min, and max of maximum onset pain: ', num2str([median(datepain(f)), min(datepain(f)), max(datepain(f))])]);
end