function ChestWallLogVDxBMIa2b
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
%fn = 'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat';
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


%for i=1:length(VDxBMICox)

v30_idx=100;
VD99BMICox = VDxBMICox(v30_idx);
cur_vd = VD99BMICox.data_exposure(:,1);
cur_vd_beta = betas(v30_idx,1);

cur_bmi = VD99BMICox.data_exposure(:,2);
cur_bmi_beta = betas(v30_idx,2);

cur_arg = [cur_vd_beta*cur_vd+cur_bmi_beta*cur_bmi];

unique_arg = unique(cur_arg);
median_arg = median(cur_arg);

rates = nan(length(unique_arg),1);
beta_vd = nan(length(unique_arg),1);
beta_bmi = nan(length(unique_arg),1);
splits = nan(length(unique_arg),1);
pvals = nan(length(unique_arg),1);
low_rates = nan(length(unique_arg),1);
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
        
        %cur_low_rate= 1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1));
        cur_low_rate = 1-sa.mSurvivalCurve{1};
        low_rates(j) = cur_low_rate(end);
        
    end
end
rates(isnan(rates)) =[];
beta_vd(isnan(beta_vd)) =[];
beta_bmi(isnan(beta_bmi)) =[];
splits(isnan(splits)) =[];
pvals(isnan(pvals)) =[];
low_rates(isnan(low_rates)) =[];
below(isnan(below)) =[];
above(isnan(above)) =[];

%% Plot volume splits for median BMI vs. low overall rates
f=figure(5);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

% volumes corresponding to splits given BMI=median(bmi)
bmi = bmi(bmi_idx);
med_bmi = median(bmi);
med_bmi_high = median(bmi(bmi>med_bmi));
med_bmi_low = median(bmi(bmi<med_bmi));


med_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi + splits./beta_vd(1);
med_high_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi_high + splits./beta_vd(1);
med_low_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi_low + splits./beta_vd(1);

h1=plot(med_bmi_low_vols,rates,'r-','LineWidth',2);hold on;
h2=plot(med_high_bmi_low_vols,rates,'g-','LineWidth',2);
h3=plot(med_low_bmi_low_vols,rates,'b-','LineWidth',2);
hold off;

%h4=plot(med_bmi_low_vols,rates,'b-','LineWidth',2);

%text(20,0.25,['V_{99,\alpha/\beta=2.1} 
ylim([0 0.3]);
xlim([0 160]);
ylabel('Lower KM overall incidence','FontSize',15);

grid on;
xlabel('Split Volume [cc]','FontSize',15);
title(['Overall CWP incidence below split vs volume'],'FontSize',15);
lgnd=legend([h3 h1 h2],...
        ['BMI = ',num2str(med_bmi_low)],...
        ['BMI = ',num2str(med_bmi)],...
        ['BMI = ',num2str(med_bmi_high)],...
        'Location','SouthEast');
 set(lgnd,'FontSize',14);

%% Plot volume splits for median BMI vs. low end rates
f=figure(6);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

% volumes corresponding to splits given BMI=median(bmi)
% bmi = bmi(bmi_idx);
% med_bmi = median(bmi);
% med_bmi_high = median(bmi(bmi>med_bmi));
% med_bmi_low = median(bmi(bmi<=med_bmi));
% 
% 
% med_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi + splits./beta_vd(1);
% med_high_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi_high + splits./beta_vd(1);
% med_low_bmi_low_vols = -(beta_bmi(1)/beta_vd(1))*med_bmi_low + splits./beta_vd(1);

h1=plot(med_bmi_low_vols,low_rates,'r-','LineWidth',2);hold on;
h2=plot(med_high_bmi_low_vols,low_rates,'g-','LineWidth',2);
h3=plot(med_low_bmi_low_vols,low_rates,'b-','LineWidth',2);

%h4=plot(med_bmi_low_vols,rates,'b-','LineWidth',2);


ylim([0 0.3]);
xlim([0 160]);
ylabel('Lower KM asymptotic incidence','FontSize',15);

grid on;
xlabel('Split Volume [cc]','FontSize',15);
title(['Asymptotic CWP incidence below split vs volume'],'FontSize',15);
lgnd=legend([h3 h1 h2],...
        ['BMI = ',num2str(med_bmi_low)],...
        ['BMI = ',num2str(med_bmi)],...
        ['BMI = ',num2str(med_bmi_high)],...
        'Location','SouthEast');
 set(lgnd,'FontSize',14);
 

 


%% plot KM splits for above rates

[rates,sort_idx] = sort(rates);
beta_vd=beta_vd(sort_idx);
beta_bmi=beta_bmi(sort_idx);
splits=splits(sort_idx);
pvals=pvals(sort_idx);


 
  % median %
[~,median_idx] = min(abs(splits-median_arg));
b_bmi=beta_bmi(median_idx);
b_vd=beta_vd(median_idx);

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
pval = sa.mpValue;
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
 text(25,0.125,['Split at median',10,'Low-risk incidence: ',...
     num2str(rates(median_idx),3),'%',10,...
     'p = ',num2str(pval,3)],...
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
 title('V_{99,a2b=2.1}+BMI at median split','fontsize',14);
 
 %% Uncertainty in parameters

f=figure(9);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

mid_range = [rates>quantile(rates,0.33)].*[rates<quantile(rates,0.66)];
mid_rates = rates(mid_range==1);
mid_pvals = pvals(mid_range==1);

sum_pval = 0;
wgtd_mid_rate = 0;
wgtd_mid_split = 0;
wgtd_mid_splits = zeros(length(mid_rates),1);

for i=1:length(mid_rates),
    [~,cur_split_idx] = min(abs(rates-mid_rates(i)));
    cur_pval = mid_pvals(i);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_mid_splits(i) = splits(cur_split_idx)*(1/cur_pval);
    wgtd_mid_split = wgtd_mid_split+splits(cur_split_idx)*(1/cur_pval);
    wgtd_mid_rate = wgtd_mid_rate+rates(cur_split_idx)*(1/cur_pval);
end
wgtd_mid_split = wgtd_mid_split/sum_pval;
wgtd_mid_rate = wgtd_mid_rate/sum_pval;
wgtd_mid_split_high =...
    wgtd_mid_split +...
    std(wgtd_mid_splits)/(sqrt(length(mid_rates))*sum_pval);

wgtd_mid_split_low =...
    wgtd_mid_split -...
    std(wgtd_mid_splits)/(sqrt(length(mid_rates))*sum_pval);



low_range = [rates<quantile(rates,0.33)];
low_rates = rates(low_range==1);
low_pvals = pvals(low_range==1);

sum_pval = 0;
wgtd_low_rate = 0;
wgtd_low_split = 0;
wgtd_low_splits = zeros(length(low_rates),1);

for i=1:length(low_rates),
    [~,cur_split_idx] = min(abs(rates-low_rates(i)));
    cur_pval = low_pvals(i);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_low_splits(i) = splits(cur_split_idx)*(1/cur_pval);
    wgtd_low_split = wgtd_low_split+splits(cur_split_idx)*(1/cur_pval);
    wgtd_low_rate = wgtd_low_rate+rates(cur_split_idx)*(1/cur_pval);
end
wgtd_low_split = wgtd_low_split/sum_pval;
wgtd_low_rate = wgtd_low_rate/sum_pval;
wgtd_low_split_high =...
    wgtd_low_split +...
    std(wgtd_low_splits)/(sqrt(length(low_rates))*sum_pval);

wgtd_low_split_low =...
    wgtd_low_split -...
    std(wgtd_low_splits)/(sqrt(length(low_rates))*sum_pval);


high_range = [rates>quantile(rates,0.66)];
high_rates = rates(high_range==1);
high_pvals = pvals(high_range==1);

sum_pval = 0;
wgtd_high_rate = 0;
wgtd_high_split = 0;
wgtd_high_splits = zeros(length(high_rates),1);

for i=1:length(high_rates),
    [~,cur_split_idx] = min(abs(rates-high_rates(i)));
    cur_pval = high_pvals(i);
    sum_pval = sum_pval+(1/cur_pval);
    wgtd_high_splits(i) = splits(cur_split_idx)*(1/cur_pval);
    wgtd_high_split = wgtd_high_split+splits(cur_split_idx)*(1/cur_pval);
    wgtd_high_rate = wgtd_high_rate+rates(cur_split_idx)*(1/cur_pval);
end
wgtd_high_split = wgtd_high_split/sum_pval;
wgtd_high_rate = wgtd_high_rate/sum_pval;
wgtd_high_split_high =...
    wgtd_high_split +...
    std(wgtd_high_splits)/(sqrt(length(high_rates))*sum_pval);

wgtd_high_split_low =...
    wgtd_high_split -...
    std(wgtd_high_splits)/(sqrt(length(high_rates))*sum_pval);


md_split=wgtd_mid_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_mid_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_mid_pct,'Color','r');
set(h_mid_pct,'LineWidth',2);


md_split=wgtd_mid_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_mid_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_mid_pct_std,'Color','r');
set(h_mid_pct_std,'LineWidth',2);
set(h_mid_pct_std,'LineStyle',':');

md_split=wgtd_mid_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_mid_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_mid_pct_std,'Color','r');
set(h_mid_pct_std,'LineStyle',':');
set(h_mid_pct_std,'LineWidth',2);

md_split=wgtd_high_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_high_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_high_pct,'Color','b');
set(h_high_pct,'LineWidth',2);

md_split=wgtd_high_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_high_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_high_pct_std,'Color','b');
set(h_high_pct_std,'LineWidth',2);
set(h_high_pct_std,'LineStyle',':');

md_split=wgtd_high_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_high_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_high_pct_std,'Color','b');
set(h_high_pct_std,'LineStyle',':');
set(h_high_pct_std,'LineWidth',2);

md_split=wgtd_low_split;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_low_pct=ezplot(bmi_vs_vd,[0,40]);
set(h_low_pct,'Color','g');
set(h_low_pct,'LineWidth',2);

md_split=wgtd_low_split_high;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_low_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_low_pct_std,'Color','g');
set(h_low_pct_std,'LineWidth',2);
set(h_low_pct_std,'LineStyle',':');

md_split=wgtd_low_split_low;
bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + md_split/b_vd;
h_low_pct_std=ezplot(bmi_vs_vd,[0,40]);
set(h_low_pct_std,'Color','g');
set(h_low_pct_std,'LineStyle',':');
set(h_low_pct_std,'LineWidth',2);


title(['Threshold for incidence in lower risk group',10,...
    'Weighted by 1/pval, errors \pm 1 S.E.'],'fontsize',15);
ylim([0 200]);
xlabel(['BMI'],'FontSize',14);
ylabel(['V_{99,a/b=2.1} [cc]'],'FontSize',14); 
grid on;

lgnd=legend([h_low_pct h_mid_pct h_high_pct],...
        ['\leq ',num2str(wgtd_low_rate*100,3),'% (Low risk)'],...        
        ['\leq ',num2str(wgtd_mid_rate*100,3),'% (Medium risk)'],...
        ['\leq ',num2str(wgtd_high_rate*100,3),'% (High risk)']);
lgnd_title=get(lgnd,'title');

set(lgnd,'fontsize',14);
set(lgnd,'location','best');


% all splits with rate < 10%
f=figure(10);  clf reset; hold on; % grid on;
set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);




end