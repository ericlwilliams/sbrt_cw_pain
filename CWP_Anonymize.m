function CWP_Anonymize
tic;
% prepare
fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
%fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
%fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');

     
     CGgrp = CGobj_current.mGrp;

    % encrypt mID 
   anon_mrns = cell(length(CGgrp),2);
   for j=1:length(CGgrp)
       
       
       tmp_id = CGgrp(j).mID;
        tmp_md5 = DataHash(tmp_id);
        anon_mrns{j,1} = tmp_id;
        anon_mrns{j,2} = tmp_md5;
       
        CGgrp(j).mID = tmp_md5;
    
   end
   
   [~,md5_inds] = sort({anon_mrns{:,2}});
   new_ids = mat2cell(md5_inds,1,ones(1,size(md5_inds,2)));
   anon_mrns(:,2) = new_ids';
   
   for k=1:length(CGgrp)
       CGgrp(k).mID = num2str(anon_mrns{k,2});
   end
   
   CGobj_current.mGrp = CGgrp;
   xlswrite(['Z:/elw/MATLAB/cw_analy/meta_data/anon_data/CWp_mrn_encrypt.xlsx'],anon_mrns)
   anon_fn=['Z:/elw/MATLAB/cw_analy/meta_data/anon_data/MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf_anon.mat'];
   
   save(anon_fn,'CGobj_current');
end
end
