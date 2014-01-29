function LoadVxDxCoxPHModelFits
tic;

%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

if isunix
    fp=strrep(fp,'G:','/media/SKI_G');
end

%fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};


CGobj = cell(length(fn),1);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

CGobj = CGobj{m};

    
    [VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
    flgCox(flganti)=false; % anti-correlations were not be considered
    VDxCox = VDxCox(flgCox);
    
    [DVxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('DVx'); % find availabe Cox models
    flgCox(flganti)=false; % anti-correlations were not be considered
    DVxCox = DVxCox(flgCox);
    
    
    vd = [VDxCox.data_exposure];
    num_vd = size(vd,2);
    vd_compdate = [VDxCox.data_hazard];
    compdate = [vd_compdate(:,1)];
    
    dv = [DVxCox.data_exposure];
    num_dv = size(dv,2); 
    
    flgcensor = [CGobj.mGrp.mFlgCensor]';
    vxdx_logl = zeros(num_vd,num_dv);
    vxdx_stats = cell(num_vd,num_dv);
    disp(['Running V_{',num2str(num_vd),...
        '} D_{',num2str(num_dv),'}']);

    for d=1:num_vd
        if ~mod(d,10)
            disp(['V_{',num2str(d),'}']);
        end
        %cur_v=1;
        for v=1:num_dv
        [~,vxdx_logl(d,v),~,vxdx_stats{d,v}]=...
            coxphfit([vd(:,d),dv(:,v)],compdate,'baseline',0,'censoring',flgcensor);            
        %cur_v=cur_v+4;
        end
    end

    %save('Z:\elw\MATLAB\cw_analy\meta_data\CW_VxDx_CoxPHM.mat','vxdx_logl','vxdx_stats');
    save('C:\Documents and Settings\williae1\cw_meta_data\RIMNER_CW_VxDx_CoxPHM.mat','vxdx_logl','vxdx_stats');
    toc;
end
