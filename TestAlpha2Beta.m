function TestAlpha2Beta

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fn=['ChestWall_a2b.mat'];

load(strcat(fp,fn));
fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1_b200.mat'];
load(strcat(fp,fn2));
CGobj=CGobj_current(1);
[VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models



screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

a2b_inf_logl = a2b_logl(1);
a2b_logl = a2b_logl(2:end);

a2b = a2b(2:end);


end