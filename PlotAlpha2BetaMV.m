function PlotAlpha2BetaMV
tic;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';


do_print = true;
if ~do_print
    disp(['NOT SAVING PLOTS'])
end
numfx = 'all';
%numfx = '5';
if isequal(numfx,'all')
    fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'];
else
    fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx',numfx,'_a2b2.1.mat'];
%    fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx',numfx,'_a2bInf.mat'];
end
%fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx3_a2b2.1.mat'];



load(strcat(fp,fn2),'CGobj_current');
CGobj = {CGobj_current};
clear CGobj_current;
CGobj = CGobj{1};

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


best_dose = 99;
%best_dose = 165;

[~,fdose_val] = min(abs(CGobj.mBinsDose - best_dose));
vd=zeros(CGobj.mNumInGrp,1); % volume v at dose d
vd(:)=0;
for k=1:CGobj.mNumInGrp
    vd(k) = CGobj.mGrp(k).fVolAtDose( CGobj.mBinsDose(fdose_val) );
end

%% V_{99}

% [cur_betas,cur_logl,~,cur_stats]=...
%     coxphfit(vd,...
%     compdate,'baseline',0,'censoring',flgcensor)
vds = vd;

vds_median = median(vds);
disp(['Median: ',num2str(vds_median,3)]);
sa_vds = repmat({classKaplanMeierCurve()},length(vds),1);

cur_i=1;
for i=1:length(vds)
    sa_vds{cur_i}.mpValue = inf;
    sa_vds{cur_i}.mHR = 0;
    
    [~,fvol] = min(abs(vds - vd(i)));
        
    flg_volbelow = vds<vds(fvol);
    if sum(flg_volbelow)< 2 || sum(flg_volbelow)>length(vds)-2,
        cur_i=cur_i+1;
        continue;
    end
    
    % Cox HR for split
    cox_beta=coxphfit(~flg_volbelow,compdate,'baseline',0,'censoring',flgcensor);
    cox_hr = exp(cox_beta);
    
    % assign properties of object sa
    survivedate={compdate(flg_volbelow); compdate(~flg_volbelow)}; % survive time of each group
    fcensor={flgcensor(flg_volbelow); flgcensor(~flg_volbelow)}; % censor flag for each group
    sa_vds{cur_i}.mSurvivalTime=survivedate;
    sa_vds{cur_i}.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa_vds{cur_i}=sa_vds{cur_i}.fCalculateSurvivalCurve();
    sa_vds{cur_i}=sa_vds{cur_i}.fCombineSurvivalTime();
    sa_vds{cur_i}=sa_vds{cur_i}.fCompareSurvivalByLogrank();
    sa_vds{cur_i}.mHR = cox_hr;
    cur_i=cur_i+1;
end

vds_pvals = [cellfun(@(x) x.mpValue, sa_vds)];

% plot pvals vs splits
fg=figure(200);  clf reset; hold on; % grid on;
set(fg,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
[sorted_vds,idx_vds] = sort(vds);hold off;
sorted_vds_pvals = vds_pvals(idx_vds);
semilogy(sorted_vds,sorted_vds_pvals);
hold on;
semilogy([min(vds) max(vds)],[0.05 0.05],'r--');
semilogy([vds_median vds_median],ylim,'g--');
semilogy([quantile(vds,0.25) quantile(vds,0.25)],ylim,'c--');
semilogy([quantile(vds,0.75) quantile(vds,0.75)],ylim,'c--');

xlabel(['V_{99} Model split value'],'FontSize',16);
ylabel(['Logrank p-value'],'FontSize',16);

% find best model (lowest logrank) in 2-3 split quartile
mid_quart_vds_pvals = sorted_vds_pvals(find([sorted_vds>quantile(sorted_vds,0.25)].*[sorted_vds<quantile(sorted_vds,0.75)]));
[min_mid_qrt_vds_pval,~] = min(mid_quart_vds_pvals);
mid_qrt_vds_best_split = find([vds_pvals==min_mid_qrt_vds_pval]);


% only print median and best
[~,vds_best_split] = min(vds_pvals);
[~,vds_median_loc] = min(abs(vds - vds_median));


%to_print = [vds_median_loc vds_best_split mid_qrt_vds_best_split] ;
to_print = [vds_median_loc];
%printing = {'Median' 'SouthEast' 'Best in 2-3rd quart'};
printing = {'Median'};


cur_sa = sa_vds(to_print);
cur_pvals = vds_pvals(to_print);
cur_splits = vds(to_print);
for j=1:length(to_print)


    
    % plot KM curves
    fg=figure(j+302);  clf reset; hold on; % grid on;
    set(fg,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    h_km(1)=stairs(cur_sa{j}.mSurvivalTimeSorted{1}./12,1-cur_sa{j}.mSurvivalCurve{1},'LineWidth',2);
    plot(cur_sa{j}.mSurvivalTimeSorted{1}(cur_sa{j}.mCensorStatistics{1}(:,1))./12,...
        1-cur_sa{j}.mSurvivalCurve{1}(cur_sa{j}.mCensorStatistics{1}(:,1)),'+','MarkerSize',12);
    h_km(2)=stairs(cur_sa{j}.mSurvivalTimeSorted{2}./12,1-cur_sa{j}.mSurvivalCurve{2},'r','LineWidth',2);
    plot(cur_sa{j}.mSurvivalTimeSorted{2}(cur_sa{j}.mCensorStatistics{2}(:,1))./12,...
        1-cur_sa{j}.mSurvivalCurve{2}(cur_sa{j}.mCensorStatistics{2}(:,1)),'r+','MarkerSize',12);
    ylim([0 0.6]);

%     str_pval = ['$p = ',num2str(cur_pvals(j),'%10.1e\n'),'$',10,...
    str_pval = ['$p \leq 0.001$',10,...
            'HR~$= ',num2str(cur_sa{j}.mHR,2),'$'];

    text(0.25,0.55,str_pval,'FontSize',20,'interpreter','latex');

    lgnd=legend(h_km,...
        strcat('$V_{99\rm{Gy}_{2.1}} <',num2str(vds_median,3),'$cc'),...
        strcat('$V_{99\rm{Gy}_{2.1}}\geq',num2str(vds_median,3),'$cc'));
    
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    set(h,'Location','NorthEast');
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Years'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    set(gca,'FontSize',18);
    
    if do_print
    set(fg,'Color','w');
    export_fig(fg,...
        [fig_loc,'km_v99_fx',numfx],'-pdf');
    end
    %title(['V_{99,a2b=2.1 GY}, Threshold: ',num2str(cur_splits(j),3), ' (',printing{j},')'],'fontsize',14);
    disp(['= V_{99, a2b=2.1}=']);
    disp(['P: ',num2str(cur_pvals(j))]);
end





%% V_{99} + BMI

[bmiCox,~,~] = CGobj.fCoxParameter_DVH('BMI');

bmi_data = [bmiCox.data_exposure];
bmi_idx = bmi_data>0;
[cur_betas,cur_logl,~,cur_stats]=...
    coxphfit([vd(bmi_idx) bmi_data(bmi_idx)],...
    compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));

disp(['= V_{99, a2b=2.1} + BMI CoxPH =']);
disp(['LogL : ',num2str(cur_logl)]);
disp(['p-vals:',10,'v_{d}: ',num2str(cur_stats.p(1)),10,'BMI: ',num2str(cur_stats.p(2))]);
disp(['Beta (V99, BMI): (',...
        num2str(cur_stats.beta(1)),',',...
        num2str(cur_stats.beta(2)),')']);
disp(['Std Err (V99, BMI): (',...
        num2str(cur_stats.se(1)),',',...
        num2str(cur_stats.se(2)),')']);

%disp(['p-val : ',num2str(cur_stats.p)]);

%KM split


g=find(bmi_idx);

vd_bmi = cur_betas(1).*vd(g) + cur_betas(2).*bmi_data(g);

vd_bmi_median = median(vd_bmi);
disp(['Median: ',num2str(vd_bmi_median)]);
sa = repmat({classKaplanMeierCurve()},length(vd_bmi),1);

% try split at each value of vd_bmi
cur_i=1;
for i=1:length(vd_bmi)
    sa{cur_i}.mpValue = inf;
    sa{cur_i}.mHR = 0;
    [~,fvol] = min(abs(vd_bmi - vd_bmi(i)));
    
    flg_volbelow = vd_bmi<vd_bmi(fvol);
    if sum(flg_volbelow)< 2 || sum(flg_volbelow)>length(vd_bmi)-2,
        cur_i=cur_i+1;
        continue;
    end
    
    % Cox HR for split
    cox_beta=coxphfit(~flg_volbelow,compdate(g),'baseline',0,'censoring',flgcensor(g));
    cox_hr = exp(cox_beta);
    
    % assign properties of object sa
    survivedate={compdate(g(flg_volbelow)); compdate(g(~flg_volbelow))}; % survive time of each group
    fcensor={flgcensor(g(flg_volbelow)); flgcensor(g(~flg_volbelow))}; % censor flag for each group
    sa{cur_i}.mSurvivalTime=survivedate;
    sa{cur_i}.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa{cur_i}=sa{cur_i}.fCalculateSurvivalCurve();
    sa{cur_i}=sa{cur_i}.fCombineSurvivalTime();
    sa{cur_i}=sa{cur_i}.fCompareSurvivalByLogrank();
    sa{cur_i}.mHR = cox_hr;
    cur_i=cur_i+1;
end

pvals = [cellfun(@(x) x.mpValue, sa)];

% plot pvals vs splits
fg=figure(100);  clf reset; hold on; % grid on;
set(fg,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
[sorted_vd_bmi,idx_vd_bmi] = sort(vd_bmi);hold off;
sorted_pvals = pvals(idx_vd_bmi);
semilogy(sorted_vd_bmi,sorted_pvals);
hold on;
semilogy([min(vd_bmi) max(vd_bmi)],[0.05 0.05],'r--');
semilogy([vd_bmi_median vd_bmi_median],ylim,'g--');
semilogy([quantile(vd_bmi,0.25) quantile(vd_bmi,0.25)],ylim,'c--');
semilogy([quantile(vd_bmi,0.75) quantile(vd_bmi,0.75)],ylim,'c--');

xlabel(['V_{99}+BMI Model split value'],'FontSize',16);
ylabel(['Logrank p-value'],'FontSize',16);

% find best model (lowest logrank) in 2-3 split quartile
mid_quart_pvals = sorted_pvals(find([sorted_vd_bmi>quantile(sorted_vd_bmi,0.25)].*[sorted_vd_bmi<quantile(sorted_vd_bmi,0.75)]));
[min_mid_qrt_pval,~] = min(mid_quart_pvals);
mid_qrt_best_split = find([pvals==min_mid_qrt_pval]);


% only print median and best
[~,best_split] = min(pvals);
[~,vd_bmi_median_loc] = min(abs(vd_bmi - vd_bmi_median));
%to_print = [vd_bmi_median_loc best_split mid_qrt_best_split] ;
%printing = {'Median' 'Best' 'Best in 2-3rd quart'};
to_print = [vd_bmi_median_loc];
printing = {'Median'};


cur_sa = sa(to_print);
cur_pvals = pvals(to_print);
cur_splits = vd_bmi(to_print);
for j=1:length(to_print)
    
    % plot KM curves
    fg=figure(j);  clf reset; hold on; % grid on;
    set(fg,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    h_km(1)=stairs(cur_sa{j}.mSurvivalTimeSorted{1}./12,1-cur_sa{j}.mSurvivalCurve{1},'LineWidth',2);
    plot(cur_sa{j}.mSurvivalTimeSorted{1}(cur_sa{j}.mCensorStatistics{1}(:,1))./12,...
        1-cur_sa{j}.mSurvivalCurve{1}(cur_sa{j}.mCensorStatistics{1}(:,1)),'+','MarkerSize',12);
    h_km(2)=stairs(cur_sa{j}.mSurvivalTimeSorted{2}./12,1-cur_sa{j}.mSurvivalCurve{2},'r','LineWidth',2);
    plot(cur_sa{j}.mSurvivalTimeSorted{2}(cur_sa{j}.mCensorStatistics{2}(:,1))./12,...
        1-cur_sa{j}.mSurvivalCurve{2}(cur_sa{j}.mCensorStatistics{2}(:,1)),'r+','MarkerSize',12);
    ylim([0 0.6]);
%     str_pval = ['$p = ',num2str(cur_pvals(j),'%10.1e\n'),'$',10,...
    str_pval = ['$p \leq 0.001$',10,...
            'HR~$= ',num2str(cur_sa{j}.mHR,2),'$'];
    %text(38,0.25,str_pval,'FontSize',16);
    text(0.25,0.55,str_pval,'FontSize',20,'interpreter','latex');
    lgnd=legend(h_km,...
        strcat('$\beta_{V_{99\rm{Gy}_{2.1}}}\times V_{99\rm{Gy}_{2.1}}+\beta_{BMI}\times BMI <',num2str(cur_splits(j),3),'$'),...
        strcat('$\beta_{V_{99\rm{Gy}_{2.1}}}\times V_{99\rm{Gy}_{2.1}}+\beta_{BMI}\times BMI \geq',num2str(cur_splits(j),3),'$'));
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    set(h,'Location','NorthEast');
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Years'],'fontsize',20);
    ylabel(['Probability of CW Pain'],'fontsize',20);
    set(gca,'FontSize',18);
    %title(['V_{99Gy,a2b=2.1 GY} + BMI, Threshold: ',num2str(cur_splits(j),3), ' (',printing{j},')'],'fontsize',14);
    if do_print
    set(fg,'Color','w');
    export_fig(fg,...
        [fig_loc,'km_v99_bmi_fx',numfx],'-pdf');
    end
end

if ~do_print
    disp(['PLOTS NOT SAVED'])
end

toc;
end