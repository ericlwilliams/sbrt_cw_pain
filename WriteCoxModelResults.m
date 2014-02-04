function WriteCoxModelResults
tic;

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b10_b200.mat',
        'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b20_b200.mat',
        'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b15_b200.mat'};
        
disp([]);

load(strcat(fp,fn{1}),'CGobj_current');
CGobj_a2b10 = CGobj_current;
load(strcat(fp,fn{2}),'CGobj_current');
CGobj_a2b20 = CGobj_current;
load(strcat(fp,fn{3}),'CGobj_current');
CGobj_a2b15 = CGobj_current;

clear CGobj_current;

%% Load DVx
[curDVxCox,flgCox,flganti] = CGobj_a2b10.fCoxParameter_DVH('DVx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b10_DVxCox = curDVxCox(flgCox);
x_a2b10_dvx=[a2b10_DVxCox.bin];

[curDVxCox,flgCox,flganti] = CGobj_a2b20.fCoxParameter_DVH('DVx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b20_DVxCox = curDVxCox(flgCox);
x_a2b20_dvx=[a2b20_DVxCox.bin];

[curDVxCox,flgCox,flganti] = CGobj_a2b15.fCoxParameter_DVH('DVx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b15_DVxCox = curDVxCox(flgCox);
x_a2b15_dvx=[a2b15_DVxCox.bin];

%% Load VDx
[curVDxCox,flgCox,flganti] = CGobj_a2b10.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b10_VDxCox = curVDxCox(flgCox);
x_a2b10_vdx=[a2b10_VDxCox.bin];

[curVDxCox,flgCox,flganti] = CGobj_a2b20.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b20_VDxCox = curVDxCox(flgCox);
x_a2b20_vdx=[a2b20_VDxCox.bin];

[curVDxCox,flgCox,flganti] = CGobj_a2b15.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
a2b15_VDxCox = curVDxCox(flgCox);
x_a2b15_vdx=[a2b15_VDxCox.bin];


fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_10-20_b200.mat'];
            
save(strcat(fp,fn),...
    'a2b10_DVxCox','a2b10_VDxCox',...
    'x_a2b10_dvx','x_a2b10_vdx',...        
    'a2b20_DVxCox','a2b20_VDxCox',...
    'x_a2b20_dvx','x_a2b20_vdx',...        
    'a2b15_DVxCox','a2b15_VDxCox',...
    'x_a2b15_dvx','x_a2b15_vdx');
%     'a2b10_DVxCox','a2b10_VDxCox',...
%     'x_a2b10_dvx','x_a2b10_vdx',...    
%     'a2b20_DVxCox','a2b20_VDxCox',...
%     'x_a2b20_dvx','x_a2b20_vdx',...        
%     'a2b15_DVxCox','a2b15_VDxCox',...
%     'x_a2b15_dvx','x_a2b15_vdx');
% ,...        
    
disp(['Saving: ',fn]);
end
