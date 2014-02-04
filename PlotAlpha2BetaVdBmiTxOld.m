function PlotAlpha2BetaVdBmiTx

tic;
screen_size=get(0,'ScreenSize');

ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
do_print =true;
fig_loc = 'Z:\elw\MATLAB\cw_analy\slides\figures\latest\';


fp = 'C:\Documents and Settings\williae1\cw_meta_data\';


fn=['ChestWall_a2b.mat'];
load(strcat(fp,fn));
%fn=['ChestWall_a2b_san.mat'];

%a2b_str = '2.1';
a2b_str = 'Inf';
fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b',a2b_str,'.mat'];

load(strcat(fp,fn2),'CGobj_current');
CGobj = {CGobj_current};
CGobj = CGobj{1};

% select patients with data
f = CGobj.fPatientsWithComplicationData();
CGobj = CGobj.fRemovePatient(~f);


screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

a2b_inf_logl = a2b_logl(1);
%% ignore inf
a2b_logl = a2b_logl(2:end);
a2b_vx = a2b_vx(2:end);
a2b = a2b(2:end);



%% Load BMI data
[bmiCox,~,~] = CGobj.fCoxParameter_DVH('BMI');
bmi_data = [bmiCox.data_exposure];
bmi_idx = bmi_data>0;

%% Load Tx data
[txCox,~,~] = CGobj.fCoxParameter_DVH('Tx');
tx_data = [txCox.data_exposure];

%% For each a2b value, find best V_{D},


% best_dose = 99;
% %best_dose = 165;
% 
% [~,fdose_val] = min(abs(CGobj.mBinsDose - best_dose));
% vd=zeros(CGobj.mNumInGrp,1); % volume v at dose d
% vd(:)=0;
% for k=1:CGobj.mNumInGrp
%     vd(k) = CGobj.mGrp(k).fVolAtDose( CGobj.mBinsDose(fdose_val) );
% end

n_grp = CGobj.mNumInGrp;
CGgrp = CGobj.mGrp;

cox_vx_logls = -inf(length(a2b),1);
cox_vx_pvals = -inf(length(a2b),1);
cox_vx_bmi_logls = -inf(length(a2b),1);
cox_vx_bmi_pvals = -inf(length(a2b),2);
cox_vx_bmi_tx_logls = -inf(length(a2b),1);
cox_vx_bmi_tx_pvals = -inf(length(a2b),3);
for i=1:length(a2b)
    % a2b_vx <- inds of best vx
    cur_best_dose = a2b_vx(i);
    cur_a2b = a2b(i);
    
    % Recalculate dose bins
    vd=zeros(n_grp,1); % volume v at dose d
    vd(:)=0;
    
    for j=1:n_grp
            
        num_fx = CGgrp(j).mFxNum;
        cur_dosebins = CGgrp(j).mDoseBins_org;
        vol = CGgrp(j).mVolCum;
        
        new_dosebins = cur_dosebins.* ...
                ((cur_a2b + ( /num_fx))/...
                (cur_a2b + 2));


        f = find( new_dosebins <= cur_best_dose ); f=f(end); % f won't be empty by definition
        if new_dosebins(f) < new_dosebins(end) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero
            vd(j) = interp1( [new_dosebins(f); new_dosebins(f+1)], [vol(f); vol(f+1)], cur_best_dose );
        else
            vd(j) = 0;
        end
    end
    
    %% Vd is current best Vd for given a2b, do multivariate
    
    %% Cox(Vd) <- sanity check
    [~,cox_vx_logls(i),~,cox_vx_stats]= coxphfit(vd,compdate,'censoring',flgcensor);
    cox_vx_pvals(i) = cox_vx_stats.p;
    
    [~,cox_vx_bmi_logls(i),~,cox_vx_bmi_stats]=coxphfit([vd(bmi_idx) bmi_data(bmi_idx)],...
    compdate(bmi_idx),'censoring',flgcensor(bmi_idx));
    cox_vx_bmi_pvals(i,1) = cox_vx_bmi_stats.p(1);
    cox_vx_bmi_pvals(i,2) = cox_vx_bmi_stats.p(2);

    [~,cox_vx_bmi_tx_logls(i),~,cox_vx_bmi_tx_stats]=coxphfit([vd(bmi_idx) bmi_data(bmi_idx) tx_data(bmi_idx)],...
    compdate(bmi_idx),'censoring',flgcensor(bmi_idx));
    cox_vx_bmi_tx_pvals(i,1) = cox_vx_bmi_tx_stats.p(1);
    cox_vx_bmi_tx_pvals(i,2) = cox_vx_bmi_tx_stats.p(2);    
    cox_vx_bmi_tx_pvals(i,3) = cox_vx_bmi_tx_stats.p(3);        

end

cox_vx_aic = 2-2.*cox_vx_logls;
cox_vx_bmi_aic = 2*2-2.*cox_vx_bmi_logls;
cox_vx_bmi_tx_aic = 2*3-2.*cox_vx_bmi_tx_logls;

cox_vx_inf_aic = 642.6;
cox_vx_bmi_inf_aic = 639.98;
cox_vx_bmi_tx_inf_aic = 638.9;


cur_fig=figure(1);  clf reset; hold on; % grid on;
set(gcf,'Position',ss_four2three);

h_vx=plot(a2b,cox_vx_aic,'r','LineWidth',2);hold on;
h_vx_bmi=plot(a2b,cox_vx_bmi_aic,'g','LineWidth',2);
h_vx_bmi_tx=plot(a2b,cox_vx_bmi_tx_aic,'b','LineWidth',2);

h_inf_vx=plot(a2b,repmat(cox_vx_inf_aic,length(a2b),1),'--r');
h_inf_vx_bmi=plot(a2b,repmat(cox_vx_bmi_inf_aic,length(a2b),1),'--g');
h_inf_vx_bmi_tx=plot(a2b,repmat(cox_vx_bmi_tx_inf_aic,length(a2b),1),'--b');

text(2,cox_vx_inf_aic+0.1,['V_{39 Gy}: ',num2str(cox_vx_inf_aic)],'Color','r','FontSize',18);
text(2,cox_vx_bmi_inf_aic+0.1,['V_{39 Gy}+BMI: ',num2str(cox_vx_bmi_inf_aic)],'Color','g','FontSize',18);
text(2,cox_vx_bmi_tx_inf_aic+0.1,['V_{39 Gy}+BMI+Tx: ',num2str(cox_vx_bmi_tx_inf_aic)],'Color','b','FontSize',18);
% h_lgnd=legend([h_vx h_inf_vx h_vx_bmi h_inf_vx_bmi h_vx_bmi_tx h_inf_vx_bmi_tx],...
%     'V_{D}','Phys: V_{39}',...
%     'V_{D}+BMI','Phys: V_{39}+BMI',...
%     'V_{D}+BMI+Tx','Phys: V_{39}+BMI+Tx');
h_lgnd=legend([h_vx h_vx_bmi h_vx_bmi_tx],...
    'V_{D}',...
    'V_{D}+BMI',...
    'V_{D}+BMI+Tx');
h_lgnd_title = get(h_lgnd,'title');
set(h_lgnd_title,'string','Cox Model','FontSize',18);
set(h_lgnd,'FontSize',16);
set(h_lgnd,'Location','NorthEast');

set(gca,'FontSize',18);
set(gca,'xminortick','on','yminortick','on');
ylabel(['AIC [2k-2ln(L)]'],'FontSize',22);
xlabel('\alpha/\beta [Gy]','fontsize',22);
%grid on;

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cox_a2b_aics'],'-pdf');
   disp(['Saving ',fig_loc,'cox_a2b_aics.pdf']);
 end
 
 %% Plot tri-variate pvalues
 
 
 
end