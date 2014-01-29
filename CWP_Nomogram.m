function CWP_Nomogram
% V99Gy_2.1 and BMI

tic;
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';
do_print=true;
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

cwp_def = 'MUTTER';

%a2b='Inf';
a2b='2.1';

%fn = {strcat(cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat')};
fn = {[cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b',a2b,'.mat']};


screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

load(strcat(fp,fn{1}),'CGobj_current');
CGobj = CGobj_current;
clear CGobj_current;




%% V30 + tx CPH

[VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
VDxCox = VDxCox(flgCox);

v30_ind=31;
if isequal(a2b,'2.1')
    v30_ind=100;
end

compdate = [VDxCox(v30_ind).data_hazard];
v30_data = [VDxCox(v30_ind).data_exposure];
flgcensor = [CGobj.mGrp.mFlgCensor]';


%% V30 + BMI CPH
[bmiCox,~,~] = CGobj.fCoxParameter_DVH('BMI');
bmi_data = [bmiCox.data_exposure];
 bmi_idx = bmi_data>0;

[cur_betas,cur_logl,cur_h,cur_stats]=...
    coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));

%% Nomogram!
[~,hzrd_2yr_ind] = min(abs(cur_h(:,1)-24));
hzrd_2yr = cur_h(hzrd_2yr_ind,2);

v30_beta = cur_betas(1);
bmi_beta = cur_betas(2);

cox_pr = (hzrd_2yr*exp(v30_beta*v30_data(bmi_idx)+bmi_beta*bmi_data(bmi_idx)));
cox_pr = cox_pr./(1+cox_pr);




end