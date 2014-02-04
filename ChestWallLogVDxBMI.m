function ChestWallLogVDxBMI
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
fn = 'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat';
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
load(strcat(fp,fn),'CGobj_current');
CGobj = CGobj_current;

% select patients with data
f = CGobj.fPatientsWithComplicationData();
CGobj = CGobj.fRemovePatient(~f);

% survival/complication time
f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(CGobj.mNumInGrp,1);
lastfollowup = inf(CGobj.mNumInGrp,1);
compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [CGobj.mGrp.mFlgCensor]';

bmi = [CGobj.mGrp.mBMI]';
bmi_idx = bmi>0;

compdate = compdate(bmi_idx);
flgcensor = flgcensor(bmi_idx);




f = cellfun(@(x) strcmpi('VDxBMI',x),CGobj.mCoxParameter(:,1)); % search the label
VDxBMICox = CGobj.mCoxParameter{f,2}; % extract Cox model result
flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), VDxBMICox); % some fields are empty or infinite, indicating no data for those values
VDxBMICox = VDxBMICox(flgCox);
tmp_beta = {VDxBMICox.beta};
is_inf = cellfun(@(x) length(x),tmp_beta);
is_inf = is_inf==1;
tmp_beta = tmp_beta(~is_inf);



VDxBMICox = VDxBMICox(~is_inf);
betas = cell2mat(tmp_beta)';

sa=classKaplanMeierCurve();

rates = nan(length(VDxBMICox),1);
beta_vd = nan(length(VDxBMICox),1);
beta_bmi = nan(length(VDxBMICox),1);
splits = nan(length(VDxBMICox),1);
pvals = nan(length(VDxBMICox),1);
flg_split2=-1;
%for i=1:length(VDxBMICox)

v30_idx=100;
cur_vd = VDxBMICox(v30_idx).data_exposure(:,1);
cur_vd_beta = betas(v30_idx,1);

cur_bmi = VDxBMICox(v30_idx).data_exposure(:,2);
cur_bmi_beta = betas(v30_idx,2);

cur_arg = [cur_vd_beta*cur_vd+cur_bmi_beta*cur_bmi];

unique_arg = unique(cur_arg);
median_arg = median(cur_arg);


rates = nan(length(unique_arg),1);
beta_vd = nan(length(unique_arg),1);
beta_bmi = nan(length(unique_arg),1);
splits = nan(length(unique_arg),1);
pvals = nan(length(unique_arg),1);
below = nan(length(unique_arg),1);
above = nan(length(unique_arg),1);

flg_split2=-1;
numstart=2;
numend=length(unique_arg)-numstart;
for j=1:length(unique_arg)
    
    
    split_arg = unique_arg(j);
    
    flg_split=cur_arg<=split_arg;
    f=length(find(flg_split));
    if f<numstart || f>numend % one group has too few patients, skip it
        continue;
    end
    if isequal(flg_split,flg_split2) % if it is the same grouping as previous, skip the computation and save the result directly
        continue;
    end
    flg_split2=flg_split;
    
    % calculate Log-Rank p-value, if sig, continue
    survivedate={compdate(flg_split); compdate(~flg_split)}; % survive time of each group
    fcensor={flgcensor(flg_split); flgcensor(~flg_split)}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    
    lr_pval = sa.mpValue;
    if lr_pval < 0.05 && cur_vd_beta>=0 && cur_bmi_beta>=0,
        
        % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
        % VDxBMImat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
        %VDxBMImat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
        
        num_below = length(survivedate{1,1});
        num_above = length(survivedate{2,1});
        cens_below = sum(fcensor{1,1});
        
        rate_below=(num_below-cens_below)/num_below;
        
        below(j) = num_below;
        above(j) = num_above;
        rates(j) = rate_below;
        beta_vd(j) = cur_vd_beta;
        beta_bmi(j) = cur_bmi_beta;
        splits(j) = split_arg;
        pvals(j) = lr_pval;
    end
end
rates(isnan(rates)) =[];
beta_vd(isnan(beta_vd)) =[];
beta_bmi(isnan(beta_bmi)) =[];
splits(isnan(splits)) =[];
pvals(isnan(pvals)) =[];
below(isnan(below)) =[];
above(isnan(above)) =[];

[rates,sort_idx] = sort(rates);
beta_vd=beta_vd(sort_idx);
beta_bmi=beta_bmi(sort_idx);
splits=splits(sort_idx);
pvals=pvals(sort_idx);

f=figure(3);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

[~,median_idx] = min(abs(splits-median_arg));
% pick best split from range around ten_pct_idx
b_bmi=beta_bmi(median_idx);b_vd=beta_vd(median_idx);md_split=splits(median_idx);
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_median_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_median_pct,'Color','k');
set(h_median_pct,'LineWidth',2);

[~,ten_pct_idx] = min(abs(rates-0.1));
% pick best split from range around ten_pct_idx
[~,best_pct_idx] = min(pvals(ten_pct_idx-10:ten_pct_idx+10));
ten_pct_idx=ten_pct_idx+best_pct_idx-1;
b_bmi=beta_bmi(ten_pct_idx);b_vd=beta_vd(ten_pct_idx);md_split=splits(ten_pct_idx);
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_ten_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_ten_pct,'Color','r');
set(h_ten_pct,'LineWidth',2);

[~,fiftn_pct_idx] = min(abs(rates-0.15));
% pick best split from range around ten_pct_idx
[~,best_pct_idx] = min(pvals(fiftn_pct_idx-10:fiftn_pct_idx+10));
fiftn_pct_idx=fiftn_pct_idx+best_pct_idx-1;
b_bmi=beta_bmi(fiftn_pct_idx);b_vd=beta_vd(fiftn_pct_idx);md_split=splits(fiftn_pct_idx);
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_fiftn_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_fiftn_pct,'Color','b');
set(h_fiftn_pct,'LineWidth',2);

[~,seven_pct_idx] = min(abs(rates-0.07));
% pick best split from range around ten_pct_idx
[~,best_pct_idx] = min(pvals(seven_pct_idx-10:seven_pct_idx+10));
seven_pct_idx=seven_pct_idx+best_pct_idx-1;
b_bmi=beta_bmi(seven_pct_idx);b_vd=beta_vd(seven_pct_idx);md_split=splits(seven_pct_idx);
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_seven_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_seven_pct,'Color','g');
set(h_seven_pct,'LineWidth',2);

title('Threshold for incidence in lower risk group','fontsize',15);
ylim([0 240]);
xlabel(['BMI'],'FontSize',14);
ylabel(['V_{30} [cc]'],'FontSize',14); 
grid on;

lgnd=legend([h_median_pct h_seven_pct h_ten_pct h_fiftn_pct],...
        ['\leq ',num2str(rates(median_idx),3),'% (Median split)'],...                
        ['\leq ',num2str(rates(seven_pct_idx),3),'% (Low risk)'],...        
        ['\leq ',num2str(rates(ten_pct_idx),3),'% (Medium risk)'],...
        ['\leq ',num2str(rates(fiftn_pct_idx),3),'% (High risk)']);
lgnd_title=get(lgnd,'title');
%set(lgnd_title,'string','Overall incidence in lower risk group:');

%set(lgnd,'interpreter','latex');
set(lgnd,'fontsize',14);
set(lgnd,'location','best');

%% plot KM splits for above rates

cur_split= splits(seven_pct_idx);
flg_split=cur_arg<=cur_split;
f=length(find(flg_split));

% calculate Log-Rank p-value, if sig, continue
survivedate={compdate(flg_split); compdate(~flg_split)}; % survive time of each group
fcensor={flgcensor(flg_split); flgcensor(~flg_split)}; % censor flag for each group
sa.mSurvivalTime=survivedate;
sa.mFlgCensor=fcensor;

sa=sa.fCalculateSurvivalCurve();
sa=sa.fCombineSurvivalTime();
sa=sa.fCompareSurvivalByLogrank();

% plot KM curves
f=figure(5);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
    1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');

 text(30,0.15,['Low risk group',10,'Low-risk incidence: ',num2str(rates(seven_pct_idx),3),'%'],...
     'FontSize',16,...
     'EdgeColor','g',...
     'LineWidth',2);
     lgnd=legend(h_km,...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \leq ',num2str(splits(seven_pct_idx),3),'$'),...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI >',num2str(splits(seven_pct_idx),3),'$'),...
        'Location','Best');

     set(lgnd,'FontSize',14);
     h=legend;
     set(h,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
  xlabel(['Months'],'fontsize',14);
  ylabel(['Probability of CW Pain'],'fontsize',14);
 title('V_{30}+BMI','fontsize',14);

 
 
% 10 %
cur_split= splits(ten_pct_idx);
flg_split=cur_arg<=cur_split;
f=length(find(flg_split));

% calculate Log-Rank p-value, if sig, continue
survivedate={compdate(flg_split); compdate(~flg_split)}; % survive time of each group
fcensor={flgcensor(flg_split); flgcensor(~flg_split)}; % censor flag for each group
sa.mSurvivalTime=survivedate;
sa.mFlgCensor=fcensor;

sa=sa.fCalculateSurvivalCurve();
sa=sa.fCombineSurvivalTime();
sa=sa.fCompareSurvivalByLogrank();

% plot KM curves
f=figure(6);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
    1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');

 text(30,0.25,['Medium risk group',10,'Low-risk incidence: ',num2str(rates(ten_pct_idx),3),'%'],...
     'FontSize',16,...
     'EdgeColor','r',...
     'LineWidth',2);
     lgnd=legend(h_km,...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \leq ',num2str(splits(ten_pct_idx),3),'$'),...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI >',num2str(splits(ten_pct_idx),3),'$'),...
        'Location','Best');

     set(lgnd,'FontSize',14);
     h=legend;
     set(h,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
  xlabel(['Months'],'fontsize',14);
  ylabel(['Probability of CW Pain'],'fontsize',14);
 title('V_{30}+BMI','fontsize',14);

 
 
 % 15 %
cur_split= splits(fiftn_pct_idx);
flg_split=cur_arg<=cur_split;
f=length(find(flg_split));

% calculate Log-Rank p-value, if sig, continue
survivedate={compdate(flg_split); compdate(~flg_split)}; % survive time of each group
fcensor={flgcensor(flg_split); flgcensor(~flg_split)}; % censor flag for each group
sa.mSurvivalTime=survivedate;
sa.mFlgCensor=fcensor;

sa=sa.fCalculateSurvivalCurve();
sa=sa.fCombineSurvivalTime();
sa=sa.fCompareSurvivalByLogrank();

% plot KM curves
f=figure(7);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
    1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');

 text(30,0.35,['High risk group',10,'Low-risk incidence: ',num2str(rates(fiftn_pct_idx),3),'%'],...
     'FontSize',16,...
     'EdgeColor','b',...
     'LineWidth',2);
     lgnd=legend(h_km,...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \leq ',num2str(splits(fiftn_pct_idx),3),'$'),...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI >',num2str(splits(fiftn_pct_idx),3),'$'),...
        'Location','Best');

     set(lgnd,'FontSize',14);
     h=legend;
     set(h,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
  xlabel(['Months'],'fontsize',14);
  ylabel(['Probability of CW Pain'],'fontsize',14);
 title('V_{30}+BMI','fontsize',14);
 
 
  % median %
cur_split= median_arg;
flg_split=cur_arg<=cur_split;
f=length(find(flg_split));

% calculate Log-Rank p-value, if sig, continue
survivedate={compdate(flg_split); compdate(~flg_split)}; % survive time of each group
fcensor={flgcensor(flg_split); flgcensor(~flg_split)}; % censor flag for each group
sa.mSurvivalTime=survivedate;
sa.mFlgCensor=fcensor;

sa=sa.fCalculateSurvivalCurve();
sa=sa.fCombineSurvivalTime();
sa=sa.fCompareSurvivalByLogrank();

% plot KM curves
f=figure(8);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
    1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
% here
%[~,median_idx]= min(abs(splits-median_arg));
 text(35,0.175,['Split at median',10,'Low-risk incidence: ',num2str(rates(median_idx),3),'%'],...
     'FontSize',16,...
     'EdgeColor','k',...
     'LineWidth',2);
     lgnd=legend(h_km,...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI \leq ',num2str(median_arg,3),'$'),...
         strcat('$\beta_{V_{30}}\times V_{30}+\beta_{BMI}\times BMI >',num2str(median_arg,3),'$'),...
        'Location','Best');

     set(lgnd,'FontSize',14);
     h=legend;
     set(h,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
  xlabel(['Months'],'fontsize',14);
  ylabel(['Probability of CW Pain'],'fontsize',14);
 title('V_{30}+BMI at median split','fontsize',14);
 
 %% Uncertainty in parameters

f=figure(9);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

sum_pval = 0;
wgtd_ten_rate = 0;
wgtd_ten_split = 0;
ten_rate_vals = [0.09:0.001:0.11];
wgtd_ten_splits = zeros(length(ten_rate_vals),1);
for i=1:length(ten_rate_vals),
    [~,cur_pct_idx] = min(abs(rates-ten_rate_vals(i)));
    cur_pval = pvals(cur_pct_idx);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_ten_splits(i) = splits(cur_pct_idx)*(1/cur_pval);
    wgtd_ten_split = wgtd_ten_split+splits(cur_pct_idx)*(1/cur_pval);
    wgtd_ten_rate = wgtd_ten_rate+rates(cur_pct_idx)*(1/cur_pval);
end
wgtd_ten_split = wgtd_ten_split/sum_pval;
wgtd_ten_rate = wgtd_ten_rate/sum_pval;
wgtd_ten_split_high =...
    wgtd_ten_split +...
    std(wgtd_ten_splits)/(sqrt(length(ten_rate_vals))*sum_pval);

wgtd_ten_split_low =...
    wgtd_ten_split -...
    std(wgtd_ten_splits)/(sqrt(length(ten_rate_vals))*sum_pval);


% seven
sum_pval = 0;
wgtd_seven_rate = 0;
wgtd_seven_split = 0;
seven_rate_vals = [0.06:0.001:0.08];
wgtd_seven_splits = zeros(length(seven_rate_vals),1);
for i=1:length(seven_rate_vals),
    [~,cur_pct_idx] = min(abs(rates-seven_rate_vals(i)));
    cur_pval = pvals(cur_pct_idx);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_seven_splits(i) = splits(cur_pct_idx)*(1/cur_pval);
    wgtd_seven_split = wgtd_seven_split+splits(cur_pct_idx)*(1/cur_pval);
    wgtd_seven_rate = wgtd_seven_rate+rates(cur_pct_idx)*(1/cur_pval);
end
wgtd_seven_split = wgtd_seven_split/sum_pval;
wgtd_seven_rate = wgtd_seven_rate/sum_pval;
wgtd_seven_split_high =...
    wgtd_seven_split +...
    std(wgtd_seven_splits)/(sqrt(length(seven_rate_vals))*sum_pval);

wgtd_seven_split_low =...
    wgtd_seven_split -...
    std(wgtd_seven_splits)/(sqrt(length(seven_rate_vals))*sum_pval);


%fifteen
sum_pval = 0;
wgtd_fiftn_split = 0;
wgtd_fiftn_rate = 0;
fiftn_rate_vals = [0.14:0.001:0.16];
wgtd_fiftn_splits = zeros(length(fiftn_rate_vals),1);
for i=1:length(fiftn_rate_vals),
    [~,cur_pct_idx] = min(abs(rates-fiftn_rate_vals(i)));
    cur_pval = pvals(cur_pct_idx);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_fiftn_splits(i) = splits(cur_pct_idx)*(1/cur_pval);
    wgtd_fiftn_split = wgtd_fiftn_split+splits(cur_pct_idx)*(1/cur_pval);
    wgtd_fiftn_rate = wgtd_fiftn_rate+rates(cur_pct_idx)*(1/cur_pval);
end
wgtd_fiftn_split = wgtd_fiftn_split/sum_pval;
wgtd_fiftn_rate = wgtd_fiftn_rate/sum_pval;
wgtd_fiftn_split_high =...
    wgtd_fiftn_split +...
    std(wgtd_fiftn_splits)/(sqrt(length(fiftn_rate_vals))*sum_pval);

wgtd_fiftn_split_low =...
    wgtd_fiftn_split -...
    std(wgtd_fiftn_splits)/(sqrt(length(fiftn_rate_vals))*sum_pval);






[~,median_idx] = min(abs(splits-median_arg));
% constant
b_bmi=beta_bmi(median_idx);
b_vd=beta_vd(median_idx);

md_split=splits(median_idx);
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_median_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_median_pct,'Color','k');
set(h_median_pct,'LineWidth',2);
set(h_median_pct,'LineStyle','--');

md_split=wgtd_ten_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_ten_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_ten_pct,'Color','r');
set(h_ten_pct,'LineWidth',2);


md_split=wgtd_ten_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_ten_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_ten_pct_std,'Color','r');
set(h_ten_pct_std,'LineWidth',2);
set(h_ten_pct_std,'LineStyle',':');

md_split=wgtd_ten_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_ten_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_ten_pct_std,'Color','r');
set(h_ten_pct_std,'LineStyle',':');
set(h_ten_pct_std,'LineWidth',2);

md_split=wgtd_fiftn_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_fiftn_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_fiftn_pct,'Color','b');
set(h_fiftn_pct,'LineWidth',2);

md_split=wgtd_fiftn_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_fiftn_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_fiftn_pct_std,'Color','b');
set(h_fiftn_pct_std,'LineWidth',2);
set(h_fiftn_pct_std,'LineStyle',':');

md_split=wgtd_fiftn_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_fiftn_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_fiftn_pct_std,'Color','b');
set(h_fiftn_pct_std,'LineStyle',':');
set(h_fiftn_pct_std,'LineWidth',2);

md_split=wgtd_seven_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_seven_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_seven_pct,'Color','g');
set(h_seven_pct,'LineWidth',2);

md_split=wgtd_seven_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_seven_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_seven_pct_std,'Color','g');
set(h_seven_pct_std,'LineWidth',2);
set(h_seven_pct_std,'LineStyle',':');

md_split=wgtd_seven_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_seven_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_seven_pct_std,'Color','g');
set(h_seven_pct_std,'LineStyle',':');
set(h_seven_pct_std,'LineWidth',2);


title(['Threshold for incidence in lower risk group',10,...
    'Weighted by 1/pval, errors \pm 1 S.E.'],'fontsize',15);
ylim([0 200]);
xlabel(['BMI'],'FontSize',14);
ylabel(['V_{30} [cc]'],'FontSize',14); 
grid on;

lgnd=legend([h_median_pct h_seven_pct h_ten_pct h_fiftn_pct],...
        ['\leq ',num2str(rates(median_idx)*100,3),'% (Median split)'],...                
        ['\leq ',num2str(wgtd_seven_rate*100,3),'% (Low risk)'],...        
        ['\leq ',num2str(wgtd_ten_rate*100,3),'% (Medium risk)'],...
        ['\leq ',num2str(wgtd_fiftn_rate*100,3),'% (High risk)']);
lgnd_title=get(lgnd,'title');

set(lgnd,'fontsize',14);
set(lgnd,'location','best');
end