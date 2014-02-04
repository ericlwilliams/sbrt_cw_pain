function CWP_IHME_etc

tic;
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';
do_print=true;
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

cwp_def = 'MUTTER';

%a2b='Inf';
a2b='2.1';

%fn = {strcat(cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat')};
fn = {[cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b',a2b,'.mat']};

%vxdx_cphm_mat_str=strcat(fp,strcat(cwp_def,'_CW_VxDx_CoxPHM.mat'));

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

load(strcat(fp,fn{1}),'CGobj_current');
CGobj = CGobj_current;
clear CGobj_current;


dpfx = [CGobj.mGrp.mDosePerFx];
tx_data  = [CGobj.mGrp.mDoseTx]';
fx = [CGobj.mGrp.mFxNum];
flgcensor = [CGobj.mGrp.mFlgCensor]';



%% V30 + tx CPH
%compdate = [txCox.data_hazard];

[VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
VDxCox = VDxCox(flgCox);

v30_ind=31;
if isequal(a2b,'2.1')
    v30_ind=100;
end

compdate = [VDxCox(v30_ind).data_hazard];
v30_data = [VDxCox(v30_ind).data_exposure];


[cur_betas,cur_logl,~,cur_stats]=...
    coxphfit([v30_data tx_data],compdate,'baseline',0,'censoring',flgcensor);

disp('V30 + Tx');
disp('(beta,se,p-val)');
disp(['v30: (',...
    num2str(cur_betas(1),3),',',...
    num2str(cur_stats.se(1)),',',...
    num2str(cur_stats.p(1))]);
disp(['TX: (',...
    num2str(cur_betas(2),3),',',...
    num2str(cur_stats.se(2)),',',...
    num2str(cur_stats.p(2))]);

v30_tx_llhd=cur_logl;

%v30_llhd = -322.65; %v30
v30_llhd = -317.87; %v30
tx_llr_test = -2*(v30_llhd)-(-2*cur_logl);
tx_llr_pval = 1-chi2cdf(tx_llr_test,1);
disp(['Tx LLRatio p-value: ',num2str(tx_llr_pval,3)]);

tx_llhd = -329.76; %v30
v30_llr_test = -2*(tx_llhd)-(-2*cur_logl);
v30_llr_pval = 1-chi2cdf(v30_llr_test,1);
disp(['V30 LLRatio p-value: ',num2str(v30_llr_pval,3)]);
disp(['LLHD: ',num2str(cur_logl,4)]);
disp(['AIC: ',num2str(2*2-2*cur_logl,5)]);



%% V30 + cm2cw CPH
[cm2cwCox,~,~] = CGobj.fCoxParameter_DVH('cm2cw');
cm2cw_data = [cm2cwCox.data_exposure];


[cur_betas,cur_logl,~,cur_stats]=...
    coxphfit([v30_data cm2cw_data],compdate,'baseline',0,'censoring',flgcensor);
disp([10]);
disp('V30 + cm2cw');
disp('(beta,se,p-val)');
disp(['v30: (',...
    num2str(cur_betas(1),3),',',...
    num2str(cur_stats.se(1)),',',...
    num2str(cur_stats.p(1))]);
disp(['cm2cw: (',...
    num2str(cur_betas(2),3),',',...
    num2str(cur_stats.se(2)),',',...
    num2str(cur_stats.p(2))]);

cm2cw_llr_test = -2*(v30_llhd)-(-2*cur_logl);
cm2cw_llr_pval = 1-chi2cdf(cm2cw_llr_test,1);
disp(['cm2cw LLRatio p-value: ',num2str(cm2cw_llr_pval,3)]);

cm2cw_llhd = -330.17; %v30
v30_llr_test = -2*(cm2cw_llhd)-(-2*cur_logl);
v30_llr_pval = 1-chi2cdf(v30_llr_test,1);
disp(['V30 LLRatio p-value: ',num2str(v30_llr_pval,3)]);
disp(['LLHD: ',num2str(cur_logl,4)]);
disp(['AIC: ',num2str(2*2-2*cur_logl,5)]);


%% V30 + BMI CPH
[bmiCox,~,~] = CGobj.fCoxParameter_DVH('BMI');
bmi_data = [bmiCox.data_exposure];
  bmi_idx = bmi_data>0;

[cur_betas,cur_logl,cur_h,cur_stats]=...
    coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
disp([10]);
disp('V30 + bmi');
disp('(beta,se,p-val)');
disp(['v30: (',...
    num2str(cur_betas(1),3),',',...
    num2str(cur_stats.se(1)),',',...
    num2str(cur_stats.p(1))]);
disp(['bmi: (',...
    num2str(cur_betas(2),3),',',...
    num2str(cur_stats.se(2)),',',...
    num2str(cur_stats.p(2))]);

v30_bmi_llhd=cur_logl;

bmi_llr_test = -2*(v30_llhd)-(-2*cur_logl);
bmi_llr_pval = 1-chi2cdf(bmi_llr_test,1);
disp(['bmi LLRatio p-value: ',num2str(bmi_llr_pval,3)]);

bmi_llhd=-335.32;
v30_llr_test = -2*(bmi_llhd)-(-2*cur_logl);
v30_llr_pval = 1-chi2cdf(v30_llr_test,1);
disp(['V30 LLRatio p-value: ',num2str(v30_llr_pval,3)]);
disp(['LLHD: ',num2str(cur_logl,4)]);
disp(['AIC: ',num2str(2*2-2*cur_logl,5)]);



return;




%% BMI + Tx
[cur_betas,cur_logl,~,cur_stats]=...
    coxphfit([bmi_data(bmi_idx) tx_data(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
disp([10]);

disp('bmi + Tx');
disp('(beta,se,p-val)');
disp(['bmi: (',...
    num2str(cur_betas(1),3),',',...
    num2str(cur_stats.se(1)),',',...
    num2str(cur_stats.p(1))]);
disp(['tx: (',...
    num2str(cur_betas(2),3),',',...
    num2str(cur_stats.se(2)),',',...
    num2str(cur_stats.p(2))]);

bmi_tx_llhd = cur_logl;%used for trivariate model

%% V30 + BMI + Tx

[cur_betas,cur_logl,~,cur_stats]=...
    coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx) tx_data(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
disp([10]);
disp(['!!!!']);
disp('V30 + bmi + Tx');
disp('(beta,se,p-val)');
disp(['v30: (',...
    num2str(cur_betas(1),3),',',...
    num2str(cur_stats.se(1)),',',...
    num2str(cur_stats.p(1))]);
disp(['bmi: (',...
    num2str(cur_betas(2),3),',',...
    num2str(cur_stats.se(2)),',',...
    num2str(cur_stats.p(2))]);
disp(['Tx: (',...
    num2str(cur_betas(3),3),',',...
    num2str(cur_stats.se(3)),',',...
    num2str(cur_stats.p(3))]);

v30_bmi_tx_llhd= cur_logl;



% v30_tx_llhd - v30_bmi_tx_llhd -> LLRD bmi
bmi_llr_test = -2*(v30_tx_llhd)-(-2*v30_bmi_tx_llhd);
bmi_llr_pval = 1-chi2cdf(bmi_llr_test,1);
disp(['BMI LLRatio p-value: ',num2str(bmi_llr_pval,3)]);

% bmi_tx_llhd - v30_bmi_tx_llhd -> LLRD v30
v30_llr_test = -2*(bmi_tx_llhd)-(-2*v30_bmi_tx_llhd);
v30_llr_pval = 1-chi2cdf(v30_llr_test,1);
disp(['V30 LLRatio p-value: ',num2str(v30_llr_pval,3)]);

% v30_bmi_llhd - v30_bmi_tx_llhd -> LLRD tx
tx_llr_test = -2*(v30_bmi_llhd)-(-2*v30_bmi_tx_llhd);
tx_llr_pval = 1-chi2cdf(tx_llr_test,1);
disp(['Tx LLRatio p-value: ',num2str(tx_llr_pval,3)]);


disp(['LLHD: ',num2str(cur_logl,4)]);
disp(['AIC: ',num2str(2*2-2*cur_logl,5)]);




% 
% cur_split = cur_betas(1)*v30_data + cur_betas(2)*cm2cw_data;
% cur_med = median(cur_split);
% 
% split_low = [cur_split<cur_med];
% 
% cox_beta=coxphfit(split_low,compdate,'baseline',0,'censoring',flgcensor);
% cox_hr = exp(cox_beta);
% 
% 

% 
% cm2cw = [CGobj.mGrp.mDistanceToChestWall]';
% 
% f = cellfun(@(x) strcmpi('CM2CW',x),CGobj.mLogRank(:,1));
% cm2cw_mat = CGobj.mLogRank{f,2};%
% 
% cm2cw_pvx = cm2cw_mat{1};
% f_pvx = cm2cw_pvx(:,6)== 1;% negative correlation (further tumor from cw -> less comps)
% unique_cm2cw=unique(cm2cw);
% unique_cm2cw=unique_cm2cw(f_pvx)
% 
% [min_pvx,idx_pvx] = min(pvx(f_pvx,5));
% 
% best_split = unique_cm2cw(idx_pvx);
% 
% sa = cm2cw_mat{2};
% sa = sa{idx_pvx};
% disp(['HR(cm2cw, min p-value) = ',num2str(sa.mHR)]);


end