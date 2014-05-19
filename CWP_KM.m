function CWP_KM
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
if isunix
    fp=strrep(fp,'G:','/media/SKI_G');
end
% fn = {'MUTTER_BC_JAN10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
%     'MUTTER_AC_JAN10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
%     'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
%km_colors={'r','b','k'};
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'};
km_colors={'k'};
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end


% survival curves
disp('Survival curves');
cur_fig=figure(1); clf reset; hold on; % grid on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

% 
%
for m=1:length(fn);
    sa = CGobj{m}.mKaplanMeierSurvivalOverall;
    stairs(sa.mSurvivalTimeSorted{1}./12,sa.mSurvivalCurve{1},km_colors{m},'LineWidth',2);
    h(m)=plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),...
       [km_colors{m},'+'],'MarkerSize',14);
end
% set(gca,'YScale','log');
set(gca,'Ylim',[0,1]);
%         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
set(gca,'xminortick','on','yminortick','on');
set(gca,'FontSize',18);
xlabel('Years','fontsize',20);
ylabel('Overall Survival','fontsize',20);
%legend(h,'Pre V_{30} (181)','Post V_{30} (135)','All (316)');
% median survival time
disp(['median survival time: ',num2str(median(sa.mSurvivalTime{1}))]);
text(3.2,0.8,['Median Survival: ',num2str(median(sa.mSurvivalTime{1})',3),'m'],'FontSize',20,'Interpreter','latex');

set(cur_fig,'Color','w');
  export_fig(cur_fig,[fig_loc,'km_os'],'-pdf');
        

f2=figure(2);clf reset;hold on;
set(f2,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);


 % assign properties of object sa
% survivedate={compdate(g(flg_dosebelow1)); compdate(g(~flg_dosebelow1))}; % survive time of each group
% fcensor={flgcensor(g(flg_dosebelow1));flgcensor(g(~flg_dosebelow1))}; % censor flag for each group
% sa.mSurvivalTime=survivedate;
% sa.mFlgCensor=fcensor;
% % compute survival curves and compare them
% sa=sa.fCalculateSurvivalCurve();
% sa=sa.fCombineSurvivalTime();
% sa=sa.fCompareSurvivalByLogrank();

sa_pre_post=classKaplanMeierCurve(); % initialize a survivalanalysis obj

sa_pre = CGobj{1}.mKaplanMeierCompOverall;
sa_post = CGobj{2}.mKaplanMeierCompOverall;

sa_pre_survivaldate=sa_pre.mSurvivalTime{1};
sa_post_survivaldate=sa_post.mSurvivalTime{1};

sa_pre_fcensor = sa_pre.mFlgCensor{1};
sa_post_fcensor = sa_post.mFlgCensor{1};

sa_pre_post_survivaldate = {[sa_pre_survivaldate];[sa_post_survivaldate]};
sa_pre_post_fcensor = {[sa_pre_fcensor];[sa_post_fcensor]};

sa_pre_post.mSurvivalTime=sa_pre_post_survivaldate;
sa_pre_post.mFlgCensor=sa_pre_post_fcensor;
sa_pre_post=sa_pre_post.fCalculateSurvivalCurve();
sa_pre_post=sa_pre_post.fCombineSurvivalTime();
sa_pre_post=sa_pre_post.fCompareSurvivalByLogrank();

pre_post_pval = sa_pre_post.mpValue;


% complication incidence curve
disp('Complication incidence curves');

for m=1:length(fn);
    sa = CGobj{m}.mKaplanMeierCompOverall;
    stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1},km_colors{m});
    g(m)=plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),...
        strcat(km_colors{m},'+'));
end
text(54,0.085,strcat('Log-rank p-value = ',num2str(pre_post_pval,3)),...
    'FontSize',12);
set(gca,'xminortick','on','yminortick','on');
set(gca,'FontSize',12);
xlabel('Months','fontsize',18);
ylabel('Probability of grade >= 2 Chestwall Pain','fontsize',18);
legend(g,'Pre V_{30} (181)','Post V_{30} (135)','All (316)','Location','SouthEast');
% median incident time
disp(['Log-rank p-value for V_{30} constraint: ',num2str(pre_post_pval,4)]);


%% Overall incidence
f3=figure(3);clf reset;hold on;
set(f3,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

sa = CGobj{3}.mKaplanMeierCompOverall;
stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},km_colors{3});
h_all=plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),...
        strcat(km_colors{3},'+'));
ylim([0 0.65]);    

%text(54,0.085,strcat('Log-rank p-value = ',num2str(pre_post_pval,3)),...
    %'FontSize',12);
set(gca,'xminortick','on','yminortick','on');
set(gca,'FontSize',16);
xlabel('Years','fontsize',20);
ylabel('Incidence of grade \geq 2 Chestwall Pain','fontsize',20);
%grid on;
%set(gca,'GridLineStyle','--')
%legend(g,'Pre V_{30} (181)','Post V_{30} (135)','All (316)','Location','SouthEast');
% median incident time
%disp(['Log-rank p-value for V_{30} constraint: ',num2str(pre_post_pval,4)]);


end
