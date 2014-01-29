function ChestWallDistanceToChestWallDisplay
tic;
% prepare
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
%fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf_b200.mat'};
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
%ss_four2three = [screen_size(3)/4 screen_size(4)/4 screen_size(3)*(3/8) screen_size(4)/2];
% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    
    CG = CGobj{m};
    
    flgcensor = [CG.mGrp.mFlgCensor]';
    cm2cw = [CG.mGrp.mDistanceToChestWall]';
    
    comp_cm2cw_avg = mean(cm2cw(~flgcensor));
    cens_cm2cw_avg = mean(cm2cw(flgcensor));
    
    [n1, xout1]=hist(cm2cw(~flgcensor),[0:1:6]);
    [n2, xout2]=hist(cm2cw(flgcensor),[0:1:6]);
    
    f(1)=figure(1);clf reset;
    set(gcf,'Position',ss_four2three);
    hold on;
    h(1)=plot(xout1,n1./sum(n1),'-ro','LineWidth',2);
    h(2)=plot([comp_cm2cw_avg comp_cm2cw_avg],ylim,'--r','LineWidth',2);
    h(3)=plot(xout2,n2./sum(n2),'-bo','LineWidth',2);
    h(4)=plot([cens_cm2cw_avg cens_cm2cw_avg],ylim,'--b','LineWidth',2);
    ylabel('Fraction of Patients','FontSize',14);
    xlabel('GTV closest distance to Chest Wall [cm]','FontSize',14);
    hold off;
   	grid on;
    h_cm2cw_lgnd = legend(h,'\geq 2 Grade CW Pain',...
                sprintf('Mean: %0.3f cc',comp_cm2cw_avg),...
                '< 2 Grade CW Pain',...
                sprintf('Mean: %0.3f cc',cens_cm2cw_avg));
    set(h_cm2cw_lgnd,'FontSize',14');
    
    %% Log-Rank p-value map
    f = cellfun(@(x) strcmpi('CM2CW',x),CGobj{m}.mLogRank(:,1));
    CM2CWmat = CG.mLogRank{f,2};%
    pvx = CM2CWmat{1};
    f = pvx(:,6)== 1;% negative correlation (further tumor from cw -> less comps)
    unique_cm2cw=unique(cm2cw);
    unique_cm2cw=unique_cm2cw(f);
    
    f(2)=figure(2); clf reset;
    set(gcf,'Position',ss_four2three);
    [min_pvx,idx_pvx] = min(pvx(f,5));
    
    [~, cm_idx] = min(abs(unique_cm2cw-1)); % index for 1cm split
    tmp_pvx = pvx(f,5);
    cm_pval = tmp_pvx(cm_idx);
    
    semilogy(unique_cm2cw,pvx(f,5),'-','LineWidth',2); hold on;
    x_cm2cw =[0:.1:max(unique_cm2cw)];
    semilogy(x_cm2cw,repmat(0.05,length(x_cm2cw),1),'r--','LineWidth',1);
    semilogy([unique_cm2cw(idx_pvx) unique_cm2cw(idx_pvx)],ylim,'g--','LineWidth',1);
    hold off; % grid on;
    xlim([0 max(unique_cm2cw)]);
    xlabel('GTV distance to Chest Wall [cm]','fontsize',14);
    ylabel('p-value','fontsize',14);
    title(strcat('Log-Rank p-values for GTV distance to chest wall'),'fontsize',14);
    str_text{1} = strcat('Min. Log-Rank p-value: ',...
            num2str(min_pvx,' %3.2e'));
    str_text{2} = strcat('at distance = ',...
                num2str(unique_cm2cw(idx_pvx)),' cm');
    
    text(2.25,0.001,str_text,'FontSize',14);

    %% KM curves for minimum p-value
    best_split = unique_cm2cw(idx_pvx);
    sa = CM2CWmat{2};
    sa = sa{idx_pvx};
    disp(['HR(cm2cw, min p-value) = ',num2str(sa.mHR)]);
    
    f(3)=figure(3);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three)
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1},'r');
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'r+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'b');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'b+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    text(40,0.3,['Log-Rank p-value: ',num2str(min_pvx,'%3.2e')],'FontSize',14);
    
    lgnd=legend(h_km,...
        strcat('Distance GTV to CW $\leq',num2str(best_split,4),'$~cm'),...
        strcat('Distance GTV to CW $>',num2str(best_split,4),'$~cm'),'Location','East');
    set(lgnd,'interpreter','latex');
    set(lgnd,'FontSize',14);

    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of CW Pain'],'fontsize',14);
    title(['GTV distance to chest wall, split at ',...
        num2str(best_split,3),' cm'],'fontsize',14);
    
    %% KM curves for split at 1cm
     %% KM curves for minimum p-value

    best_split = unique_cm2cw(cm_idx);
    sa = CM2CWmat{2};
    sa = sa{cm_idx};
    disp(['HR(cm2cw, 1cm) = ',num2str(sa.mHR)]);
    
    f(4)=figure(4);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three)
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'r');
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'r+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'b');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'b+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    
    ylim([0 0.65]);
    textbp(['Log-Rank p-value: ',num2str(cm_pval,'%3.2e'),10,...
        'HR = ',num2str(sa.mHR,3)],'FontSize',16);
    
    lgnd=legend(h_km,...
        strcat('Distance GTV to CW $\leq 1$~cm'),...
        strcat('Distance GTV to CW $> 1$~cm'),'Location','Best');
    set(lgnd,'interpreter','latex');
    set(lgnd,'FontSize',18);
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'fontsize',18);
    xlabel(['Years'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    title(['GTV distance to chest wall, split at 1 cm'],'fontsize',20);
    
end