function WriteCoxModelResults2

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_300-500_b200.mat'];% 300, 400, 500
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_75-200_b200.mat'];% 75, 100, 200
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_30-50_b200.mat'];
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_10-20_b200.mat'];
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_3-5_b200.mat'];
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_Inf-2_b200.mat'];
load(strcat(fp,fn))
fn=['MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_b200.mat'];
save(strcat(fp,fn))

end