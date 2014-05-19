function CWP_ROC
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

flgcensor = [CGobj{m}.mGrp.mFlgCensor]';
flgcomp = ~flgcensor;
for m = 1:length(fn)
    
    V30=zeros(CGobj{m}.mNumInGrp,1);
    
    for k=1:CGobj{m}.mNumInGrp
        V30(k) = CGobj{m}.mGrp(k).fVolAtDose( CGobj{m}.mBinsDose(30));
    end
    
   v30_roc([V30 flgcomp]);     
end
end