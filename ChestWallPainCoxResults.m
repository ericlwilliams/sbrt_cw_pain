function ChestWallPainCoxResults
tic;
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';
do_print=true;
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fig_num=1;
if isunix
    fp=strrep(fp,'G:','/media/SKI_G');
end
%cwp_def = 'RIMNER';
cwp_def = 'MUTTER';

fn = {strcat(cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat')};
%fn = {strcat(cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat')};

vxdx_cphm_mat_str=strcat(fp,strcat(cwp_def,'_CW_VxDx_CoxPHM.mat'));
%vxdx_cphm_mat_str='C:\Documents and Settings\williae1\cw_meta_data\MUTTER_CW_VxDx_CoxPHM.mat';


%'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
    %'RIMNER_BC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
    %'RIMNER_AC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
    
    %'MUTTER_BC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat',...
    %'MUTTER_AC_SEPT10_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};

CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
ss_full = screen_size;
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end


for m = 1:length(fn)
    disp(fn{m});
    uni_print={};
    
    mv_logl = ones(8);% Tx, Dose/Fx, cm2cw, NumFx, BMI, vd
    mv_aic = zeros(8);
    mv_pvals = ones(8,8,2);% Tx, Dose/Fx, cm2cw, NumFx, BMI, vd
    
    % Tx
    [txCox,~,~] = CGobj{m}.fCoxParameter_DVH('Tx');
    tx_uni_p = txCox.p;
    tx_uni_logl=txCox.logl;
    cur_print= {'Tx' txCox.beta txCox.se...
        txCox.logl txCox.p};
    uni_print = [uni_print;cur_print];
    mv_logl(1,1) = txCox.logl;
    mv_aic(1,1) = -2*txCox.logl + 2;
    mv_pvals(1,1,1) = tx_uni_p;
    
    
    % Dose per fraction
    [dosePerFxCox,~,~] = CGobj{m}.fCoxParameter_DVH('DperFx');
    dperfx_uni_p=dosePerFxCox.p;
    dperfx_uni_logl=dosePerFxCox.logl;
    cur_print = {'Dose/Fx' dosePerFxCox.beta dosePerFxCox.se...
        dosePerFxCox.logl dosePerFxCox.p};
    uni_print = [uni_print;cur_print];
    mv_logl(2,2) = dosePerFxCox.logl;
    mv_aic(2,2) = -2*dosePerFxCox.logl + 2;
    mv_pvals(2,2,1) = dperfx_uni_p;
    
     % cm2cw
    [cm2cwCox,~,~] = CGobj{m}.fCoxParameter_DVH('cm2cw');
    cm2cw_uni_p = cm2cwCox.p;
    cm2cw_uni_logl=cm2cwCox.logl;
    cur_print = {'cm2cw' cm2cwCox.beta cm2cwCox.se...
        cm2cwCox.logl cm2cwCox.p};
    uni_print = [uni_print;cur_print];
    mv_logl(3,3) = cm2cwCox.logl;
    mv_aic(3,3) = -2*cm2cwCox.logl + 2;
    mv_pvals(3,3,1) = cm2cw_uni_p;
    
    % NumFx
    [numFxCox,~,~] = CGobj{m}.fCoxParameter_DVH('Fx');
    numfx_uni_p=numFxCox.p;
    numfx_uni_logl=numFxCox.logl;
    cur_print = {'NumFx' numFxCox.beta numFxCox.se...
        numFxCox.logl numFxCox.p};
   uni_print = [uni_print;cur_print];
    mv_logl(4,4) = numFxCox.logl;
    mv_aic(4,4) = -2*numFxCox.logl + 2;
    mv_pvals(4,4,1) = numfx_uni_p;
    
   % BMI
    [bmiCox,~,~] = CGobj{m}.fCoxParameter_DVH('BMI');
    bmi_uni_p=bmiCox.p;
    bmi_uni_logl=bmiCox.logl;
    cur_print = {'BMI' bmiCox.beta bmiCox.se...
        bmiCox.logl bmiCox.p};
    uni_print = [uni_print;cur_print];
    mv_logl(5,5) = bmiCox.logl;
    mv_aic(5,5) = -2*bmiCox.logl + 2;
    mv_pvals(5,5,1) = bmi_uni_p;
   
    % KPS
    [curCox,~,~] = CGobj{m}.fCoxParameter_DVH('KPS');
    cur_print = {'KPS' curCox.beta curCox.se...
        curCox.logl curCox.p};
    uni_print = [uni_print;cur_print];
  
     % Sex
    [curCox,~,~] = CGobj{m}.fCoxParameter_DVH('Gender');
    cur_print = {'Sex' curCox.beta curCox.se...
        curCox.logl curCox.p};
    uni_print = [uni_print;cur_print];
    
    % Age
    [curCox,~,~] = CGobj{m}.fCoxParameter_DVH('Age');
    cur_print = {'AgeAtTx' curCox.beta curCox.se...
        curCox.logl curCox.p};
    uni_print = [uni_print;cur_print];
     
      

    
    
    %% D_{Vx} Cox PH Model results
    [DVxCox,flgCox,flganti] = CGobj{m}.fCoxParameter_DVH('DVx'); % find availabe Cox models
    flgCox(flganti)=false; % anti-correlations were not be considered
    DVxCox = DVxCox(flgCox);

   
        
    logl = [DVxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
    [mx,doseloc]=max(logl); % the best fitting of Cox model
    lowCI68 = mx - 0.5; % 68% confidence
    lowCI95 = mx - 1.96; % 95% confidence
    
    % from meta
    %      lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
    %        lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
    
    % mutter best
    cur_print = {'D_{83}' DVxCox(84).beta DVxCox(84).se...
        DVxCox(84).logl DVxCox(84).p};
    uni_print = [uni_print;cur_print];       
    
    % Rimner best
    cur_print = {'D_{16}' DVxCox(17).beta DVxCox(17).se...
        DVxCox(17).logl DVxCox(17).p};
    uni_print = [uni_print;cur_print];       
    
    
    d83_data = [DVxCox(84).data_exposure];
    d16_data = [DVxCox(17).data_exposure];
        
    p = ones(size(DVxCox));
    p = [DVxCox.p]';
    [min_p,ploc] = min(p);
    
    mv_logl(7,7) = mx;
    mv_aic(7,7) = -2*mx + 2;
    mv_pvals(7,7,1)=min_p;

     x_dvx=CGobj{m}.mBinsVol(flgCox);    
          
     disp(['Best D_{V} Cox Model at D_{',num2str(x_dvx(ploc)),'} is:']);
     disp(DVxCox(ploc));
%     
%   
    
%here dv
    % logl
    cur_fig=figure(fig_num); clf reset; 
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    plot(x_dvx, [DVxCox.logl],'x-','LineWidth',2,'MarkerSize',10);
    fig_num=fig_num+1;
    %set(gcf,'Position',ss_four2three);
    hold on; 
    ci68=plot(x_dvx,repmat(lowCI68,size(x_dvx)),'r--'); 
    ci95=plot(x_dvx,repmat(lowCI95,size(x_dvx)),'c--');
    h_mx_logl=plot([x_dvx(doseloc) x_dvx(doseloc)],ylim,'g--','LineWidth',2);
    xlim([0 max(x_dvx)]);
    set(gca,'fontsize',18);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    
    xlabel('(D_{V}) Volume [cc]','fontsize',24);
    ylabel('CPH log-likelihood','fontsize',24);
    logl_str = {['Max LogL = ',10, num2str(mx,5),' at D_{',num2str(x_dvx(doseloc),3),' cc}']};
    lgnd=legend([ci68,ci95, h_mx_logl],['Low 68% CI','Low 95% CI',logl_str],'Location','Best');

    set(lgnd,'FontSize',20);
    set(lgnd,'Location','NorthEast');
    if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_dv_llhds'],'-pdf');
    disp(['Saving ',fig_loc,'cph_dv_llhds.pdf']);
    end

    
    %% DVx Correlatinos
    %dvx_corrs = corr([DVxCox.data_exposure],[DVxCox.data_exposure]).^2;
    dvx_corrs = corr([DVxCox.data_exposure],[DVxCox.data_exposure]);
    cur_fig=figure(1000);
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);    
    imagesc(1:length(dvx_corrs),1:length(dvx_corrs),dvx_corrs);
    set(gca,'YDir','normal');
    colorbar;
    title('D_{V} R Correlations','FontSize',18);
    ylabel('(D_{V}) Volume [cc]','FontSize',24)
    xlabel('(D_{V}) Volume [cc]','FontSize',24)

    if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_vd_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'cph_vd_dv_r.pdf']);
    end

    
    % p-vals
    
    cur_fig=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);    
    semilogy(x_dvx,[DVxCox.p],'x-','LineWidth',2,'MarkerSize',10); 
    hold on; 
    semilogy(x_dvx,repmat(0.05,size(x_dvx)),'r--'); 
    h_min_pval = semilogy([x_dvx(ploc) x_dvx(ploc)],ylim,'g--','LineWidth',2);
    hold off;
    xlim([0 max(x_dvx)]);
    set(gca,'fontsize',18);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    pval_str = {['Min p-val = ', num2str(min_p,2),' at D_{',num2str(x_dvx(ploc),4),' cc}']};
    lgnd=legend(h_min_pval,pval_str,'Location','Best');
    set(lgnd,'FontSize',20);
    xlabel('(D_{V}) Volume [cc]','fontsize',24);
    ylabel('CPH p-value','fontsize',24);
    if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_dv_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'cph_dv_pvals.pdf']);
    end

    disp(['best fit by log likelihood and p-value: ',num2str([doseloc,ploc])]);
    
    %             doseloc = 31;
    %DVxCox = DVxCox(doseloc);
    disp('Best Cox Model:');
    disp(DVxCox(doseloc));
    dv_doseloc=doseloc;
    
    %% V_{Dx}
        
    [VDxCox,flgCox,flganti] = CGobj{m}.fCoxParameter_DVH('VDx'); % find availabe Cox models
    flgCox(flganti)=false; % anti-correlations were not be considered
    VDxCox = VDxCox(flgCox);
    logl = [VDxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
    [mx,doseloc]=max(logl); % the best fitting of Cox model
    lowCI68 = mx - 0.5; % 68% confidence
    lowCI95 = mx - 1.96; % 95% confidence
    
    
    
    
     %% Vdx Correlatinos
    vdx_corrs = corr([VDxCox.data_exposure],[VDxCox.data_exposure]).^2;
    figure(1001);
    imagesc(1:length(vdx_corrs),1:length(vdx_corrs),vdx_corrs);
    set(gca,'YDir','normal');
    colorbar;
    title('V_{D} R^2 Correlations','FontSize',14);
    xlabel('(V_{D}) Dose [Gy]','FontSize',14)
    ylabel('(V_{D}) Dose [Gy]','FontSize',14)
       
   
    
    
   
     %% Vdx Correlatinos
    dvx_vdx_corrs = corr([DVxCox.data_exposure],[VDxCox.data_exposure]).^2;
    cur_fig=figure(100001);
    %imagesc(1:length(dvx_vdx_corrs),1:length(dvx_vdx_corrs),dvx_vdx_corrs);
    imagesc(dvx_vdx_corrs);
    set(gca,'YDir','normal');
    colorbar;
    title('R(V_{D},D_{V}) Correlations','FontSize',20);
    ylabel('(D_{V}) Volume [cc]','FontSize',20)
    xlabel('(V_{D}) Dose [Gy]','FontSize',20)
    set(gca,'FontSize',18);
    
          if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_vd_dv_r'],'-pdf');
    disp(['Saving ',fig_loc,'cph_vd_dv_r.pdf']);
          end

    
    v42_data = [VDxCox(43).data_exposure];
    
    v39_data = [VDxCox(40).data_exposure];
    v30_data = [VDxCox(31).data_exposure];
    
    %% only for a2b = 2.1 and a2b=0
    
    v99_data = [VDxCox(40).data_exposure];
    v165_data = [VDxCox(40).data_exposure];
    %v99_data = [VDxCox(100).data_exposure];
    %v165_data = [VDxCox(166).data_exposure];
    
    %mutter best
    cur_print = {'V_{39}' VDxCox(40).beta VDxCox(40).se...
        VDxCox(40).logl VDxCox(40).p};
    uni_print = [uni_print;cur_print]; 
    
    % Rimner best
    cur_print = {'V_{42}' VDxCox(43).beta VDxCox(43).se...
        VDxCox(43).logl VDxCox(43).p};
    uni_print = [uni_print;cur_print]; 
    
    cur_print = {'V_{30}' VDxCox(31).beta VDxCox(31).se...
        VDxCox(31).logl VDxCox(31).p};
    uni_print = [uni_print;cur_print]; 
    
    cur_print = {['V_{' num2str(doseloc-1) '}'] VDxCox(doseloc).beta VDxCox(doseloc).se...
        VDxCox(doseloc).logl VDxCox(doseloc).p};
    uni_print = [uni_print;cur_print]; 
    
    cur_print = {['D_{' num2str(dv_doseloc-1) '}'] DVxCox(dv_doseloc).beta DVxCox(dv_doseloc).se...
        DVxCox(dv_doseloc).logl DVxCox(dv_doseloc).p};
    
    uni_print = [uni_print;cur_print]; 
    
    
    
    
    %             num = cellfun(@(x) size(x,1),{VDxCox.data_exposure});
    x_dvx=CGobj{m}.mBinsDose(flgCox);    
    
    % logl
    cur_fig=figure(fig_num); clf reset; 
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    plot(x_dvx, [VDxCox.logl],'x-','LineWidth',2,'MarkerSize',10);
    fig_num=fig_num+1;
    %set(gcf,'Position',ss_four2three);
    hold on; 
    ci68=plot(x_dvx,repmat(lowCI68,size(x_dvx)),'r--'); 
    ci95=plot(x_dvx,repmat(lowCI95,size(x_dvx)),'c--');
    h_mx_logl=plot([x_dvx(doseloc) x_dvx(doseloc)],ylim,'g--','LineWidth',2);
    xlim([0 max(x_dvx)]);
    set(gca,'fontsize',18);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    
    xlabel('(V_{D}) Dose [Gy]','fontsize',24);
    ylabel('CPH log-likelihood','fontsize',24);
    logl_str = {['Max LogL = ',10, num2str(mx,5),' at D_{',num2str(x_dvx(doseloc),4),' Gy}']};
    lgnd=legend([ci68,ci95, h_mx_logl],['Low 68% CI','Low 95% CI',logl_str],'Location','Best');
    set(lgnd,'FontSize',20);

    if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_vd_llhds'],'-pdf');
    disp(['Saving ',fig_loc,'cph_vd_llhds.pdf']);
    end

    % p-vals
    p = ones(size(VDxCox));
    p(flgCox) = [VDxCox.p]';
    [min_p,ploc] = min(p);
    
    
    
    
    cur_fig=figure(fig_num); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    fig_num=fig_num+1;
    %set(gcf,'Position',ss_four2three);
    semilogy(x_dvx,[VDxCox.p],'x-','LineWidth',2,'MarkerSize',10); 
    hold on; 
    semilogy(x_dvx,repmat(0.05,size(x_dvx)),'r--'); 
    h_min_pval = semilogy([x_dvx(ploc) x_dvx(ploc)],ylim,'g--','LineWidth',2);
    hold off;
    xlim([0 max(x_dvx)]);
    set(gca,'fontsize',18);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    pval_str = {['Min p-val = ', num2str(min_p,2),' at V_{',num2str(x_dvx(ploc),4),' Gy}']};
    lgnd=legend(h_min_pval,pval_str,'Location','Best');
    set(lgnd,'FontSize',20);
    xlabel('(V_{D}) Dose [Gy]','fontsize',24);
    ylabel('CPH p-value','fontsize',24);
    
    if do_print
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'cph_vd_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'cph_vd_pvals.pdf']);
    end

    disp(['best fit by log likelihood and p-value: ',num2str([x_dvx(doseloc),x_dvx(ploc)])]);
    
    %             doseloc = 31;
    disp('Best Cox Model:');
    disp(VDxCox(doseloc));
            
    
    
    %% Multi-variate models
    


    compdate = [txCox.data_hazard];
    compdate = [compdate(:,1)];
    flgcensor = [CGobj{m}.mGrp.mFlgCensor]';
    
    % Tx + Dose/Fx
    tx_data = [txCox.data_exposure];
    dosePerFx_data = [dosePerFxCox.data_exposure];
    [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data dosePerFx_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(1,2) = cur_logl;
    mv_aic(1,2) = -2*cur_logl + 2*2;
    mv_pvals(1,2,:) = cur_stats.p;
    
    
    % Tx + cm2cw
    cm2cw_data = [cm2cwCox.data_exposure];
    [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data cm2cw_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(1,3) = cur_logl; 
    mv_aic(1,3) = -2*cur_logl + 2*2;
    mv_pvals(1,3,:) = cur_stats.p;
    
    % Tx + numFx
    numFx_data = [numFxCox.data_exposure];
    [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data numFx_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(1,4) = cur_logl; 
    mv_aic(1,4) = -2*cur_logl + 2*2;
    mv_pvals(1,4,:) = cur_stats.p;
    
    % Tx + BMI
    bmi_data = [bmiCox.data_exposure];
    bmi_idx = bmi_data>0;
    [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) bmi_data(bmi_idx)],...
            compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    mv_logl(1,5) = cur_logl; 
    mv_aic(1,5) = -2*cur_logl + 2*2;
    mv_pvals(1,5,:) = cur_stats.p;
    
    % V42 + KPS
    [kpsCox,~,~] = CGobj{m}.fCoxParameter_DVH('KPS');
    kps_data = [kpsCox.data_exposure];
    [~,cur_logl,~,cur_stats]=...
            coxphfit([kps_data tx_data],...
            compdate,'baseline',0,'censoring',flgcensor);
    
    
    
    
    % dosePerFx + cm2cw
    [~,cur_logl,~,cur_stats]=...
            coxphfit([dosePerFx_data cm2cw_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(2,3) = cur_logl; 
    mv_aic(2,3) = -2*cur_logl + 2*2;
    mv_pvals(2,3,:) = cur_stats.p;
    
    % dosePerFx + numFx
    [~,cur_logl,~,cur_stats]=...
            coxphfit([dosePerFx_data numFx_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(2,4) = cur_logl; 
    mv_aic(2,4) = -2*cur_logl + 2*2;
    mv_pvals(2,4,:) = cur_stats.p;
    
    % dosePerFx + bmi
    [~,cur_logl,~,cur_stats]=...
            coxphfit([dosePerFx_data(bmi_idx) bmi_data(bmi_idx)],...
            compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    mv_logl(2,5) = cur_logl; 
    mv_aic(2,5) = -2*cur_logl + 2*2;
    mv_pvals(2,5,:) = cur_stats.p;
    
    % cm2cw + numFx
    [~,cur_logl,~,cur_stats]=...
            coxphfit([cm2cw_data numFx_data],compdate,'baseline',0,'censoring',flgcensor);
    mv_logl(3,4) = cur_logl; 
    mv_aic(3,4) = -2*cur_logl + 2*2;
    mv_pvals(3,4,:) = cur_stats.p;
    
    % cm2cw + bmi
    [~,cur_logl,~,cur_stats]=...
            coxphfit([cm2cw_data(bmi_idx) bmi_data(bmi_idx)],...
            compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    mv_logl(3,5) = cur_logl; 
    mv_aic(3,5) = -2*cur_logl + 2*2;
    mv_pvals(3,5,:) = cur_stats.p;
    
    % numFx + bmi
    [~,cur_logl,~,cur_stats]=...
            coxphfit([numFx_data(bmi_idx) bmi_data(bmi_idx)],...
            compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    mv_logl(4,5) = cur_logl; 
    mv_aic(4,5) = -2*cur_logl + 2*2;
    mv_pvals(4,5,:) = cur_stats.p;
    
   
    
    %% V_{39} + D_{83}
    % D_{83} models stored in DVxCox
    % V_{39} models stored in VDxCox
    % Run bivariate Cox PH Model with V_D and D_V
      
  
    % Run LoadVxDxCoxPHModelFits.m first then 
    % loads vxdx_logl and vxdx_stats
    load(vxdx_cphm_mat_str)
    
    figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    % shift x and y axis -1 to account for mBinsDose/mBinsVol shift
    y_vxdx = [0:size(vxdx_logl',1)-1];
    x_vxdx = [0:size(vxdx_logl',2)-1];
    imagesc(x_vxdx,y_vxdx,vxdx_logl');
    set(gca,'YDir','normal');
    colorbar;
    title('V_{D} + D_{V} Cox PH Model Log-Likelihood','FontSize',14);
    xlabel('(V_{D}) Dose [Gy]','FontSize',14)
    ylabel('(D_{V}) Volume [cc]','FontSize',14)
   
    [vxdx_max_dv, vxdx_max_vd] =find(vxdx_logl'==max(max(vxdx_logl')));
    vxdx_max_logl=max(max(vxdx_logl'));
    mv_logl(6,7)=vxdx_max_logl;
    mv_aic(6,7) = -2*vxdx_max_logl + 2*2;

    
    % p-values
    vd_pvals = cellfun(@(x) x.p(1), vxdx_stats);
    dv_pvals = cellfun(@(x) x.p(2), vxdx_stats);
    
    vd_vxdx_min_pval = vd_pvals(vxdx_max_vd);
    dv_vxdx_min_pval = dv_pvals(vxdx_max_dv);
   
    mv_pvals(6,7,2)=vd_vxdx_min_pval;
    mv_pvals(6,7,1)=dv_vxdx_min_pval;
    
    idx_v30=31;
    mv_pvals(7,8,2) = min(vd_pvals(idx_v30,:));
    mv_pvals(7,8,1) = min(dv_pvals(idx_v30,:));
    mv_logl(7,8) =  max(vxdx_logl(30,:));% V_{30)}

    disp(['']);
    disp('*******');
     disp(['Max LogL for D_{83} + V_{39} Model at: D_{',...
        num2str(vxdx_max_dv),...
        '} + V_{',...
        num2str(vxdx_max_vd),...
        '}']);
    disp(['Corresponding pvalues, D_{',...
        num2str(vxdx_max_dv),...
        '}: ',...
        num2str(dv_vxdx_min_pval),...
        ' and V_{',...
        num2str(vxdx_max_vd),...
        '}: ',...
        num2str(vd_vxdx_min_pval)]);
    
    % P-values for V_{39} coefficient in V_{39} + D_{83} model
    figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    imagesc(x_vxdx,y_vxdx,vd_pvals');
    set(gca,'YDir','normal');
    colorbar
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    ylabel(['(D_{V}) Volume [Gy]'],'FontSize',14);
    xlabel(['(V_{D}) Dose [Gy]'],'FontSize',14);
    title('P-values for V_{D} coefficient in V_{D} + D_{V} CPHM','FontSize',14);
    
    
    % P-values for D_{83} coefficient in V_{39} + D_{83} model
    figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    imagesc(x_vxdx,y_vxdx,dv_pvals');
    set(gca,'YDir','normal');
    colorbar
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    ylabel(['(D_{V}) Volume [Gy]'],'FontSize',14);
    xlabel(['(V_{D}) Dose [Gy]'],'FontSize',14);
    title('P-values for D_{V} coefficient in V_{D} + D_{V} CPHM','FontSize',14);
    
    
    
    % Maximum p-val between the two
    use_vx_pvals = vd_pvals>dv_pvals;
    max_vxdx_pvals = vd_pvals.*use_vx_pvals + dv_pvals.*~use_vx_pvals;
    
    figure(fig_num);
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    h_vxdx_pvals = imagesc(x_vxdx,y_vxdx,max_vxdx_pvals');
    set(gca,'YDir','normal');
    colorbar;
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    %h_vxdx_pvals_axis = get(h_vxdx_pvals,'Parent');
    %clim_vxdx_pvals = get(h_vxdx_pvals_axis,'CLim');
    ylabel(['(D_{V}) Volume [cc]'],'FontSize',14);
    xlabel(['(V_{D}) Dose [Gy]'],'FontSize',14);
    title(['Max P-value between \beta_{\rm{V}_{\rm{D}}} and \beta_{\rm{D}_{\rm{V}}}'],'FontSize',14);
    
    
    % Average when both are less than 0.05
    
    avg_vxdx_pvals = vd_pvals+dv_pvals;
    avg_vxdx_pvals = avg_vxdx_pvals./2;
    [min_avg_row,min_avg_col] = find(avg_vxdx_pvals==min(min(avg_vxdx_pvals)));

    disp(['Minimum Avg D_{83}/V_{39} combination at V_{',...
    num2str(min_avg_row),...
    '} (p <= ',...
        num2str(vd_pvals(min_avg_row,min_avg_col)),...
        ') D_{',...
        num2str(min_avg_col),...
        '} (p <= ',...
        num2str(dv_pvals(min_avg_row,min_avg_col)),...
        ')']);
    
    use_vx = vd_pvals<=0.05;
    use_dx = dv_pvals<=0.05;
    use_min = use_vx.*use_dx;
    min_pvals = vd_pvals.*use_min + dv_pvals.*use_min;
    min_pvals = min_pvals./2;
    min_pvals(~min_pvals)=Inf;
    
    figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    h_vxdx_min_avg_pvals = imagesc(x_vxdx,y_vxdx,min_pvals');
    set(gca,'YDir','normal');
    colorbar;
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    %h_vxdx_pvals_axis = get(h_vxdx_pvals,'Parent');
    %clim_vxdx_pvals = get(h_vxdx_pvals_axis,'CLim');
    ylabel(['(D_{V}) Volume [cc]'],'FontSize',14);
    xlabel(['(V_{D}) Dose [Gy] '],'FontSize',14);
    title(['Mean p-value for both < 0.05'],'FontSize',14);
    


    %% Number of Fractions and V_{39}
    
    [biCox,flgBiCox,flgBiAnti] = CGobj{m}.fCoxParameter_DVH('VDxFx');
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('VDx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    dosebins = CGobj{m}.mBinsDose(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
    mv_logl(6,6)=max(uni_logl);
    mv_aic(6,6) = -2*max_bi_logl + 2;
    mv_logl(4,6)=max_bi_logl;
    mv_aic(4,6) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(dosebins,[uni_logl' repmat(numfx_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni V_{39}','Uni NumFx','V_{39} + NumFx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at V_{',...
        num2str(dosebins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('V_{D} + NumFx Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    vd_vd_numfx_min_pval = bi_p(1,idx_max_bi_logl);
    numfx_vd_numfx_min_pval = bi_p(2,idx_max_bi_logl);
    
    
    mv_pvals(8,8,1) = uni_p(idx_v30);% V_{30)}
    mv_logl(8,8) =  uni_logl(idx_v30);% V_{30)}
    
    mv_pvals(4,8,1) = bi_p(2,idx_v30);
    mv_pvals(4,8,2) = bi_p(1,idx_v30);
    mv_logl(4,8) =  bi_logl(idx_v30);% V_{30)}
    
    [mv_pvals(6,6,1),idx_min_uni_p]  = min(uni_p);
    uniCox = uniCox(flgUniBiCox);
    disp(['Best V_{39} Cox Model at V_{',num2str(dosebins(idx_min_uni_p)),'} is:']);
    disp(uniCox(idx_min_uni_p));
    
    mv_pvals(4,6,1)  = numfx_vd_numfx_min_pval;
    mv_pvals(4,6,2)  = vd_vd_numfx_min_pval;
            
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(dosebins,[uni_p' bi_p' repmat(numfx_uni_p,size(uni_p'))],'.-');
    hold on; semilogy(dosebins,repmat(0.05,size(dosebins)),'g--','LineWidth',3); 
    h5=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);set(h1(4),'MarkerSize',12);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni V_{39}','Bi V_{39}','Uni NumFx','Bi NumFx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
    title('V_{D} + Number of Fractions');
    
    
     %% V_{39} + Tx
        
    %[biCox,flgBiCox,flgBiAnti] = CGobj{m}.fCoxParameter_DVH('VDxTx');
    f = cellfun(@(x) strcmpi('VDxTx',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('VDx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    dosebins = CGobj{m}.mBinsDose(flgUniBiCox);
    
     % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
    mv_logl(1,6)=max_bi_logl;
    mv_aic(1,6) = -2*max_bi_logl + 2*2;
    
    
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(dosebins,[uni_logl' repmat(tx_uni_logl,size(uni_logl')) bi_logl'],'.-');
    h_mx_bi_logl=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni V_{39}','Uni Tx','V_{39} + Tx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at V_{',...
        num2str(dosebins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('V_{D} + Tx Uni- and Bi-variate Log-likelihood');
    

    %plot p-values
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd and tx p-values at maximum logl point
    vd_vd_tx_min_pval = bi_p(1,idx_max_bi_logl);
    tx_vd_tx_min_pval = bi_p(2,idx_max_bi_logl);
    
    mv_pvals(1,6,1)  = tx_vd_tx_min_pval;
    mv_pvals(1,6,2)  = vd_vd_tx_min_pval;
    
    mv_pvals(1,8,1) = bi_p(2,idx_v30);
    mv_pvals(1,8,2) = bi_p(1,idx_v30);
    mv_logl(1,8) =  bi_logl(idx_v30);% V_{30)}
    

    f3=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h3=semilogy(dosebins,[uni_p' bi_p' repmat(tx_uni_p,size(uni_p'))],'.-');
    hold on; semilogy(dosebins,repmat(0.05,size(dosebins)),'g--','LineWidth',3); 
    h5=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
        hold off;
    xlim([0 max(dosebins)]);
    set(h3(1),'Color','b');set(h3(1),'LineWidth',2);set(h3(1),'MarkerSize',12)
    set(h3(2),'Color','b');set(h3(2),'LineWidth',2);set(h3(2),'LineStyle',':');set(h3(2),'MarkerSize',12)
    set(h3(3),'Color','r');set(h3(3),'LineWidth',2);set(h3(3),'MarkerSize',12);set(h3(3),'LineStyle',':');
    set(h3(4),'Color','r');set(h3(4),'LineWidth',2);set(h3(4),'MarkerSize',12)
    legend([h3(1) h3(2) h3(4) h3(3) h5],'Uni V_{39}','Bi V_{39}','Uni Tx','Bi Tx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');

    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
    title('V_{D} + Prescription Dose');
 
    %% V_{39} + cm2cw

    f = cellfun(@(x) strcmpi('VDxCM2CW',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('VDx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    dosebins = CGobj{m}.mBinsDose(flgUniBiCox);
    
     % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
    mv_logl(3,6)=max_bi_logl;
    mv_aic(3,6) = -2*max_bi_logl + 2*2;
    
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(dosebins,[uni_logl' repmat(cm2cw_uni_logl,size(uni_logl')) bi_logl'],'.-');
    h_mx_bi_logl=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni V_{39}','Uni cm2cw','V_{39} + cm2cw','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at V_{',...
        num2str(dosebins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('V_{D} + Distance to CW Uni- and Bi-variate Log-likelihood');
    

    %plot p-values
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd and tx p-values at maximum logl point
    vd_vd_cm2cw_min_pval = bi_p(1,idx_max_bi_logl);
    cm2cw_vd_cm2cw_min_pval = bi_p(2,idx_max_bi_logl);
    
    mv_pvals(3,6,1)  = cm2cw_vd_cm2cw_min_pval;
    mv_pvals(3,6,2)  = vd_vd_cm2cw_min_pval;
    
    mv_pvals(3,8,1) = bi_p(2,idx_v30);
    mv_pvals(3,8,2) = bi_p(1,idx_v30);
    mv_logl(3,8) =  bi_logl(idx_v30);% V_{30)}
    
    f3=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h3=semilogy(dosebins,[uni_p' bi_p' repmat(cm2cw_uni_p,size(uni_p'))],'.-');
    hold on; semilogy(dosebins,repmat(0.05,size(dosebins)),'g--','LineWidth',3); 
    h5=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
        hold off;
    xlim([0 max(dosebins)]);
    set(h3(1),'Color','b');set(h3(1),'LineWidth',2);set(h3(1),'MarkerSize',12)
    set(h3(2),'Color','b');set(h3(2),'LineWidth',2);set(h3(2),'LineStyle',':');set(h3(2),'MarkerSize',12)
    set(h3(3),'Color','r');set(h3(3),'LineWidth',2);set(h3(3),'MarkerSize',12);set(h3(3),'LineStyle',':');
    set(h3(4),'Color','r');set(h3(4),'LineWidth',2);set(h3(4),'MarkerSize',12)
    legend([h3(1) h3(2) h3(4) h3(3) h5],'Uni V_{39}','Bi V_{39}','Uni cm2cw','Bi cm2cw','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');

    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
    title('V_{D} + Distance to CW');
     %% V_{39} + BMI
        
    %[biCox,flgBiCox,flgBiAnti] = CGobj{m}.fCoxParameter_DVH('VDxBMI');
    f = cellfun(@(x) strcmpi('VDxBMI',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('VDx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox; % cut off when univariate goes anti-correlated
    
    dosebins = CGobj{m}.mBinsDose(flgUniBiCox);
    
     % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);

    mv_logl(5,6)=max_bi_logl;
    mv_aic(5,6) = -2*max_bi_logl + 2*2;

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(dosebins,[uni_logl' repmat(bmi_uni_logl,size(uni_logl')) bi_logl'],'.-');
    h_mx_bi_logl=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni V_{39}','Uni BMI','V_{39} + BMI','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at V_{',...
        num2str(dosebins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('V_{D} + BMI Uni- and Bi-variate Log-likelihood');
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
     % vd and tx p-values at maximum logl point
    vd_vd_bmi_min_pval = bi_p(1,idx_max_bi_logl);
    bmi_vd_bmi_min_pval = bi_p(2,idx_max_bi_logl);
    
    mv_pvals(5,6,1)  = bmi_vd_bmi_min_pval;
    mv_pvals(5,6,2)  = vd_vd_bmi_min_pval;
    
    mv_pvals(5,8,1) = bi_p(2,idx_v30);
    mv_pvals(5,8,2) = bi_p(1,idx_v30);
    mv_logl(5,8) =  bi_logl(idx_v30);% V_{30)}
    
    f4=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h4=semilogy(dosebins,[uni_p' bi_p' repmat(bmi_uni_p,size(uni_p'))],'.-');
    hold on; semilogy(dosebins,repmat(0.05,size(dosebins)),'g--','LineWidth',3); 
    h5=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    
    hold off;
    xlim([0 max(dosebins)]);
    set(h4(1),'Color','b');set(h4(1),'LineWidth',2);set(h4(1),'MarkerSize',12)
    set(h4(2),'Color','b');set(h4(2),'LineWidth',2);set(h4(2),'LineStyle',':');set(h4(2),'MarkerSize',12)
    set(h4(3),'Color','r');set(h4(3),'LineWidth',2);set(h4(3),'MarkerSize',12);set(h4(3),'LineStyle',':');
    set(h4(4),'Color','r');set(h4(4),'LineWidth',2);set(h4(4),'MarkerSize',12)
    
    legend([h4(1) h4(2) h4(4) h4(3) h5],'Uni V_{39}','Bi V_{39}','Uni BMI','Bi BMI','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
    title('V_{D} + BMI');
    
     %% V_{39} + Dose/fraction
        
    %[biCox,flgBiCox,flgBiAnti] = CGobj{m}.fCoxParameter_DVH('VDxBMI');
    f = cellfun(@(x) strcmpi('VDxDperFx',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('VDx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox; % cut off when univariate goes anti-correlated
    
    dosebins = CGobj{m}.mBinsDose(flgUniBiCox);
    
     % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);

    mv_logl(2,6)=max_bi_logl;
    mv_aic(2,6) = -2*max_bi_logl + 2*2;

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(dosebins,[uni_logl' repmat(dperfx_uni_logl,size(uni_logl')) bi_logl'],'.-');
    h_mx_bi_logl=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(dosebins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni V_{39}','Uni Dose/Fx','V_{39} + Dose/Fx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at V_{',...
        num2str(dosebins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('V_{D} + Dose/Fx Uni- and Bi-variate Log-likelihood');
        
    
    % plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
        
    vd_vd_dperfx_min_pval = bi_p(1,idx_max_bi_logl);
    dperfx_vd_dperfx_min_pval = bi_p(2,idx_max_bi_logl);
    
    mv_pvals(2,6,1)  = dperfx_vd_dperfx_min_pval;
    mv_pvals(2,6,2)  = vd_vd_dperfx_min_pval;
    
    mv_pvals(2,8,1) = bi_p(2,idx_v30);
    mv_pvals(2,8,2) = bi_p(1,idx_v30);
    mv_logl(2,8) =  bi_logl(idx_v30);% V_{30)}
    
    f5=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h5=semilogy(dosebins,[uni_p' bi_p' repmat(dperfx_uni_p,size(uni_p'))],'.-');
    hold on; semilogy(dosebins,repmat(0.05,size(dosebins)),'g--','LineWidth',3); 
    h6=plot([dosebins(idx_max_bi_logl) dosebins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    hold off;
    xlim([0 max(dosebins)]);
    set(h5(1),'Color','b');set(h5(1),'LineWidth',2);set(h5(1),'MarkerSize',12)
    set(h5(2),'Color','b');set(h5(2),'LineWidth',2);set(h5(2),'LineStyle',':');set(h5(2),'MarkerSize',12)
    set(h5(3),'Color','r');set(h5(3),'LineWidth',2);set(h5(3),'MarkerSize',12);set(h5(3),'LineStyle',':');
    set(h5(4),'Color','r');set(h5(4),'LineWidth',2);set(h5(4),'MarkerSize',12)
    
    legend([h5(1) h5(2) h5(4) h5(3) h6],'Uni V_{39}','Bi V_{39}','Uni Dose/Fx','Bi Dose/Fx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    xlabel('(V_{D}) Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
    title('V_{D} + Dose per Fraction');
   
  %% DV + Number of Fractions
  
    [biCox,flgBiCox,flgBiAnti] = CGobj{m}.fCoxParameter_DVH('DVxFx');
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
    
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('DVx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    volbins = CGobj{m}.mBinsVol(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
     mv_logl(4,7)=max_bi_logl;
     mv_aic(4,7) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(volbins,[uni_logl' repmat(numfx_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni D_{83}','Uni NumFx','D_{83} + NumFx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at D_{',...
        num2str(volbins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(D_{83}) Volume [cc]','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('D_{83} + NumFx Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    dv_dv_numfx_min_pval = bi_p(1,idx_max_bi_logl);
    numfx_dv_numfx_min_pval = bi_p(2,idx_max_bi_logl);
    
     mv_pvals(4,7,1)  = numfx_dv_numfx_min_pval;
     mv_pvals(4,7,2)  = dv_dv_numfx_min_pval;
    
        
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(volbins,[uni_p' bi_p' repmat(numfx_uni_p,size(uni_p'))]);
    hold on; semilogy(volbins,repmat(0.05,size(volbins)),'g--','LineWidth',3); 
    h5=plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);%set(h1(1),'MarkerSize',5)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');%set(h1(2),'MarkerSize',5)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'LineStyle',':');%set(h1(3),'MarkerSize',5);
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);%set(h1(4),'MarkerSize',5);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni D_{83}','Bi D_{83}','Uni NumFx','Bi NumFx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);
    title('D_{V} + Number of Fractions');
    
    
     %% DV + Prescription Dose
  
    f = cellfun(@(x) strcmpi('DVxTx',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
        
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('DVx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    volbins = CGobj{m}.mBinsVol(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
     mv_logl(1,7)=max_bi_logl;
     mv_aic(1,7) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(volbins,[uni_logl' repmat(tx_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni D_{83}','Uni Tx','D_{83} + Tx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at D_{',...
        num2str(volbins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('D_{V} + Tx Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    dv_dv_tx_min_pval = bi_p(1,idx_max_bi_logl);
    tx_dv_tx_min_pval = bi_p(2,idx_max_bi_logl);
    
     mv_pvals(1,7,1)  = tx_dv_tx_min_pval;
     mv_pvals(1,7,2)  = dv_dv_tx_min_pval;
        
        
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(volbins,[uni_p' bi_p' repmat(tx_uni_p,size(uni_p'))]);
    hold on; semilogy(volbins,repmat(0.05,size(volbins)),'g--','LineWidth',3); 
    h5=plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);%set(h1(1),'MarkerSize',5)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');%set(h1(2),'MarkerSize',5)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'LineStyle',':');%set(h1(3),'MarkerSize',5);
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);%set(h1(4),'MarkerSize',5);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni D_{V}','Bi D_{V}','Uni Tx','Bi Tx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);
    title('D_{V} + Prescription Dose');
    
      %% DV + Distance to Chest Wall
  
    f = cellfun(@(x) strcmpi('DVxCM2CW',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
        
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('DVx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    volbins = CGobj{m}.mBinsVol(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
     mv_logl(3,7)=max_bi_logl;
     mv_aic(3,7) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(volbins,[uni_logl' repmat(cm2cw_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni D_{V}','Uni cm2cw','D_{V} + cm2cw','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at D_{',...
        num2str(volbins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('D_{V} + Distance to CW Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    dv_dv_cm2cw_min_pval = bi_p(1,idx_max_bi_logl);
    cm2cw_dv_cm2cw_min_pval = bi_p(2,idx_max_bi_logl);
    
     mv_pvals(3,7,1)  = cm2cw_dv_cm2cw_min_pval;
     mv_pvals(3,7,2)  = dv_dv_cm2cw_min_pval;
        
        
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(volbins,[uni_p' bi_p' repmat(cm2cw_uni_p,size(uni_p'))]);
    hold on; semilogy(volbins,repmat(0.05,size(volbins)),'g--','LineWidth',3); 
    h5=plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);%set(h1(1),'MarkerSize',5)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');%set(h1(2),'MarkerSize',5)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'LineStyle',':');%set(h1(3),'MarkerSize',5);
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);%set(h1(4),'MarkerSize',5);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni D_{V}','Bi D_{V}','Uni cm2cw','Bi cm2cw','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);
    title('D_{V} + cm2cw');
    
        %% DV + BMI
  
    f = cellfun(@(x) strcmpi('DVxBMI',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
        
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('DVx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    volbins = CGobj{m}.mBinsVol(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
     mv_logl(5,7)=max_bi_logl;
     mv_aic(5,7) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(volbins,[uni_logl' repmat(bmi_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni D_{V}','Uni BMI','D_{V} + BMI','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at D_{',...
        num2str(volbins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('D_{V} + Distance to CW Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    dv_dv_bmi_min_pval = bi_p(1,idx_max_bi_logl);
    bmi_dv_bmi_min_pval = bi_p(2,idx_max_bi_logl);
    
     mv_pvals(5,7,1)  = bmi_dv_bmi_min_pval;
     mv_pvals(5,7,2)  = dv_dv_bmi_min_pval;

     
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(volbins,[uni_p' bi_p' repmat(bmi_uni_p,size(uni_p'))]);
    hold on; semilogy(volbins,repmat(0.05,size(volbins)),'g--','LineWidth',3); 
    h5=plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);%set(h1(1),'MarkerSize',5)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');%set(h1(2),'MarkerSize',5)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'LineStyle',':');%set(h1(3),'MarkerSize',5);
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);%set(h1(4),'MarkerSize',5);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni D_{V}','Bi D_{V}','Uni BMI','Bi BMI','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);
    title('D_{V} + BMI');
    
    %% DV + Dose/Fraction
  
    f = cellfun(@(x) strcmpi('DVxDperFx',x),CGobj{m}.mCoxParameter(:,1)); % search the label
    biCox = CGobj{m}.mCoxParameter{f,2}; % extract Cox model result
    flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), biCox); % some fields are empty or infinite, indicating no data for those values
    tmp_beta = {biCox.beta};
    is_inf = cellfun(@(x) length(x),tmp_beta);
    is_inf = is_inf==1;
    tmp_beta = tmp_beta(~is_inf);
    f = cell2mat(tmp_beta)';
    flgBiAnti = f<0;
    
    flgBiCox(flgBiAnti(:,1))=false; % anti-correlations were not be considered
        
    [uniCox,flgUniCox,flgUniAnti] = CGobj{m}.fCoxParameter_DVH('DVx');
    flgUniCox(flgUniAnti)=false; % anti-correlations were not be considered
    
    %flgUniBiCox=flgUniCox&flgBiCox;
    flgUniBiCox=flgUniCox;
    
    volbins = CGobj{m}.mBinsVol(flgUniBiCox);
    
    % plot log-likelihood
    uni_logl = [uniCox(flgUniBiCox).logl];
    bi_logl = [biCox(flgUniBiCox).logl];
    
    [max_bi_logl,idx_max_bi_logl] = max(bi_logl);
    
     mv_logl(2,7)=max_bi_logl;
     mv_aic(2,7) = -2*max_bi_logl + 2*2;
    

    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    
    set(gcf,'Position',ss_four2three);
    hold on;
    h1=plot(volbins,[uni_logl' repmat(dperfx_uni_logl,size(uni_logl')) bi_logl'],'.-');
   plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);set(h1(1),'MarkerSize',12)
    set(h1(2),'Color','r');set(h1(2),'LineWidth',2);set(h1(2),'MarkerSize',12)
    set(h1(3),'Color','k');set(h1(3),'LineWidth',2);set(h1(3),'MarkerSize',12);set(h1(3),'LineStyle',':');
    
    legend([h1(1) h1(2) h1(3)],'Uni D_{V}','Uni Dose/Fx','D_{V} + Dose/Fx','Location','Best');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    bi_logl_str = {['Max LogL = ', num2str(max_bi_logl,5),' at D_{',...
        num2str(volbins(idx_max_bi_logl),4)],'}'};
    text(2,-337,bi_logl_str,'fontsize',16);
    hold off;
    grid on;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('Log-likelihood','fontsize',18);
    title('D_{V} + Distance to CW Uni- and Bi-variate Log-likelihood');
    
    
    %plot pvals
    bi_p = [biCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    uni_p = [uniCox(flgUniBiCox).p]; %ignore anti-correlated for V_{Dx}
    
    % vd adn numfx p-values at maximum logl point
    dv_dv_dperfx_min_pval = bi_p(1,idx_max_bi_logl);
    dperfx_dv_dperfx_min_pval = bi_p(2,idx_max_bi_logl);
    
     mv_pvals(2,7,1)  = dperfx_dv_dperfx_min_pval;
     mv_pvals(2,7,2)  = dv_dv_dperfx_min_pval;
        
        
    f1=figure(fig_num); clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h1=semilogy(volbins,[uni_p' bi_p' repmat(dperfx_uni_p,size(uni_p'))]);
    hold on; semilogy(volbins,repmat(0.05,size(volbins)),'g--','LineWidth',3); 
    h5=plot([volbins(idx_max_bi_logl) volbins(idx_max_bi_logl)],...
        ylim,'m--','LineWidth',2);
    xlim([0 max(volbins)]);
    set(h1(1),'Color','b');set(h1(1),'LineWidth',2);%set(h1(1),'MarkerSize',5)
    set(h1(2),'Color','b');set(h1(2),'LineWidth',2);set(h1(2),'LineStyle',':');%set(h1(2),'MarkerSize',5)
    set(h1(3),'Color','r');set(h1(3),'LineWidth',2);set(h1(3),'LineStyle',':');%set(h1(3),'MarkerSize',5);
    set(h1(4),'Color','r');set(h1(4),'LineWidth',2);%set(h1(4),'MarkerSize',5);
    
    legend([h1(1) h1(2) h1(4) h1(3) h5],'Uni D_{V}','Bi D_{V}','Uni Dose/Fx','Bi Dose/Fx','Max LogL','Location','SouthWest');
    set(gca,'fontsize',14);
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    hold off;
    xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);
    title('D_{V} + Dose/Fx');
    
    
    
    %% print pvals
    
    print_pvals = mv_pvals;
    
    % re-arrange

    %V_{39} 6->1
    print_pvals(1,1,:) = mv_pvals(6,6,:);
    print_pvals(1,2,:) = mv_pvals(6,8,:);% VD + V30
    print_pvals(1,3,:) = mv_pvals(6,7,:);% VD + DV
    print_pvals(1,4,:) = mv_pvals(1,6,:);% VD + Tx
    print_pvals(1,5,:) = mv_pvals(2,6,:);% VD + DperFx
    print_pvals(1,6,:) = mv_pvals(3,6,:);% VD + cm2cw
    print_pvals(1,7,:) = mv_pvals(4,6,:);% VD + NumFx
    print_pvals(1,8,:) = mv_pvals(5,6,:);% VD + BMI
    
    %V_{30} 8->2
    print_pvals(2,2,:) = mv_pvals(8,8,:);
    print_pvals(2,3,:) = mv_pvals(7,8,:);% D83 + V30
    print_pvals(2,4,:) = mv_pvals(1,8,:);% Tx + V30
    print_pvals(2,5,:) = mv_pvals(2,8,:);% VD + DperFx
    print_pvals(2,6,:) = mv_pvals(3,8,:);% VD + cm2cw
    print_pvals(2,7,:) = mv_pvals(4,8,:);% VD + NumFx
    print_pvals(2,8,:) = mv_pvals(5,8,:);% VD + BMI
    
    %D83 7->3
    print_pvals(3,3,:) = mv_pvals(7,7,:);% D83 + V30
    print_pvals(3,4,:) = mv_pvals(1,7,:);% Tx + V30
    print_pvals(3,5,:) = mv_pvals(2,7,:);% VD + DperFx
    print_pvals(3,6,:) = mv_pvals(3,7,:);% VD + cm2cw
    print_pvals(3,7,:) = mv_pvals(4,7,:);% VD + NumFx
    print_pvals(3,8,:) = mv_pvals(5,7,:);% VD + BMI
    
    % Tx 1->4
    print_pvals(4,4,:) = mv_pvals(1,1,:);% D83 + V30
    print_pvals(4,5,:) = mv_pvals(1,2,:);% Tx + V30
    print_pvals(4,6,:) = mv_pvals(1,3,:);% VD + DperFx
    print_pvals(4,7,:) = mv_pvals(1,4,:);% VD + cm2cw
    print_pvals(4,8,:) = mv_pvals(1,5,:);% VD + NumFx
    
    % D/Fx 2->5
    print_pvals(5,5,:) = mv_pvals(2,2,:);% D83 + V30
    print_pvals(5,6,:) = mv_pvals(2,3,:);% Tx + V30
    print_pvals(5,7,:) = mv_pvals(2,4,:);% VD + DperFx
    print_pvals(5,8,:) = mv_pvals(2,5,:);% VD + cm2cw
    
    % cm2cw 3->6
    print_pvals(6,6,:) = mv_pvals(3,3,:);% D83 + V30
    print_pvals(6,7,:) = mv_pvals(3,4,:);% Tx + V30
    print_pvals(6,8,:) = mv_pvals(3,5,:);% VD + DperFx

    % numFx 4->7
    print_pvals(7,7,:) = mv_pvals(4,4,:);% D83 + V30
    print_pvals(7,8,:) = mv_pvals(4,5,:);% Tx + V30
    
    % BMI 5->8
    print_pvals(8,8,:) = mv_pvals(5,5,:);% D83 + V30
    
        
    mv_pvals=print_pvals;
    
    mv_logl(mv_logl==1)=NaN;
    mv_aic(mv_aic==0)=NaN;
    mv_pvals(mv_pvals==1)=NaN;
   
    
    figure(fig_num);clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_two2two);
    set(gcf,'Position',ss_full);
    avg_pvals = (mv_pvals(:,:,1)+mv_pvals(:,:,2))./2;
    h_pvals = imagesc(avg_pvals);
    set(gca,'YDir','normal');
    set(gca,'XTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    set(gca,'YTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    set(gca,'FontSize',18);
    set(h_pvals,'alphadata',~isnan(avg_pvals));
    colorbar;
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    
    
    first_pvals = mv_pvals(:,:,1);
    diag_pvals = diag(first_pvals);
    first_pvals = first_pvals-diag(diag(first_pvals));
    first_pvals(first_pvals==0)=NaN;
    first_isnan = isnan(first_pvals);
    %first_pvals(first_pvals<0.01)=0.01;
    first_pvals = first_pvals(~first_isnan);
    
    second_pvals = mv_pvals(:,:,2);
    second_isnan = isnan(second_pvals);
    %second_pvals(second_pvals<0.01)=0.01;
    second_pvals=second_pvals(~second_isnan);
    
    
    first_pvals_strs = num2str(first_pvals(:),'%2.2g');
    second_pvals_strs = num2str(second_pvals(:),'%2.2g');
    diag_pvals_strs = num2str(diag_pvals(:),'%2.2g');
    
    first_pvals_strs = strtrim(cellstr(first_pvals_strs));%# Remove any space padding
    second_pvals_strs = strtrim(cellstr(second_pvals_strs));%# Remove any space padding
    diag_pvals_strs = strtrim(cellstr(diag_pvals_strs));%# Remove any space padding
    
    [x_pvals,y_pvals] = meshgrid(1:8); %# x/y coordiantes for strings
    x_diag_pvals = [1:8];
    y_diag_pvals = [1:8];
    
    
    text(x_pvals(~first_isnan),y_pvals(~first_isnan)+.25,first_pvals_strs,...
            'HorizontalAlignment','right','FontSize',16);
    text(x_pvals(~second_isnan),y_pvals(~second_isnan)-.25,second_pvals_strs,...
            'HorizontalAlignment','left','FontSize',16);
    text(x_diag_pvals,y_diag_pvals,diag_pvals_strs,...
            'HorizontalAlignment','center','FontSize',20);

    % diagnal lines to split boxes
    line([1 8.5],[0 7.5],'Color','k');
    line([2 8.5],[0 6.5],'Color','k');
    line([3 8.5],[0 5.5],'Color','k');
    line([4 8.5],[0 4.5],'Color','k');
    line([5 8.5],[0 3.5],'Color','k');
    line([6 8.5],[0 2.5],'Color','k');
    line([7 8.5],[0 1.5],'Color','k');
    
    % vertical lines
    line([0.01 0.01],ylim,'Color','k');
    line(xlim,[0.01 0.01],'Color','k');
    line([1.5 1.5],[0 2.5],'Color','k');
    line([2.5 2.5],[0 3.5],'Color','k');
    line([3.5 3.5],[0 4.5],'Color','k');    
    line([4.5 4.5],[0 5.5],'Color','k');    
    line([5.5 5.5],[0 6.5],'Color','k');    
    line([6.5 6.5],[0 7.5],'Color','k');    
    line([7.5 7.5],[0 8.5],'Color','k');    
    
    % horizontal lines
    line(xlim,[0 0],'Color','k');
    line([0.5 5.5],[1.5 1.5],'Color','k');
    line([1.5 5.5],[2.5 2.5],'Color','k');
    line([2.5 5.5],[3.5 3.5],'Color','k');
    line([3.5 5.5],[4.5 4.5],'Color','k');
    line([4.5 5.5],[5.5 5.5],'Color','k');
    line([5.5 6.5],[6.5 6.5],'Color','k');
    line([6.5 7.5],[7.5 7.5],'Color','k');
    line([7.5 8.5],[8.5 8.5],'Color','k');
    title('Multi-var Cox PH Model p-values','FontSize',18);
    
    
      print_llhd = mv_logl;
    
    % re-arrange

    %V_{39} 6->1
    print_llhd(1,1,:) = mv_logl(6,6,:);
    print_llhd(1,2,:) = mv_logl(6,8,:);% VD + V30
    print_llhd(1,3,:) = mv_logl(6,7,:);% VD + DV
    print_llhd(1,4,:) = mv_logl(1,6,:);% VD + Tx
    print_llhd(1,5,:) = mv_logl(2,6,:);% VD + DperFx
    print_llhd(1,6,:) = mv_logl(3,6,:);% VD + cm2cw
    print_llhd(1,7,:) = mv_logl(4,6,:);% VD + NumFx
    print_llhd(1,8,:) = mv_logl(5,6,:);% VD + BMI
    
    %V_{30} 8->2
    print_llhd(2,2,:) = mv_logl(8,8,:);
    print_llhd(2,3,:) = mv_logl(7,8,:);% D83 + V30
    print_llhd(2,4,:) = mv_logl(1,8,:);% Tx + V30
    print_llhd(2,5,:) = mv_logl(2,8,:);% VD + DperFx
    print_llhd(2,6,:) = mv_logl(3,8,:);% VD + cm2cw
    print_llhd(2,7,:) = mv_logl(4,8,:);% VD + NumFx
    print_llhd(2,8,:) = mv_logl(5,8,:);% VD + BMI
    
    %D83 7->3
    print_llhd(3,3,:) = mv_logl(7,7,:);% D83 + V30
    print_llhd(3,4,:) = mv_logl(1,7,:);% Tx + V30
    print_llhd(3,5,:) = mv_logl(2,7,:);% VD + DperFx
    print_llhd(3,6,:) = mv_logl(3,7,:);% VD + cm2cw
    print_llhd(3,7,:) = mv_logl(4,7,:);% VD + NumFx
    print_llhd(3,8,:) = mv_logl(5,7,:);% VD + BMI
    
    % Tx 1->4
    print_llhd(4,4,:) = mv_logl(1,1,:);% D83 + V30
    print_llhd(4,5,:) = mv_logl(1,2,:);% Tx + V30
    print_llhd(4,6,:) = mv_logl(1,3,:);% VD + DperFx
    print_llhd(4,7,:) = mv_logl(1,4,:);% VD + cm2cw
    print_llhd(4,8,:) = mv_logl(1,5,:);% VD + NumFx
    
    % D/Fx 2->5
    print_llhd(5,5,:) = mv_logl(2,2,:);% D83 + V30
    print_llhd(5,6,:) = mv_logl(2,3,:);% Tx + V30
    print_llhd(5,7,:) = mv_logl(2,4,:);% VD + DperFx
    print_llhd(5,8,:) = mv_logl(2,5,:);% VD + cm2cw
    
    % cm2cw 3->6
    print_llhd(6,6,:) = mv_logl(3,3,:);% D83 + V30
    print_llhd(6,7,:) = mv_logl(3,4,:);% Tx + V30
    print_llhd(6,8,:) = mv_logl(3,5,:);% VD + DperFx

    % numFx 4->7
    print_llhd(7,7,:) = mv_logl(4,4,:);% D83 + V30
    print_llhd(7,8,:) = mv_logl(4,5,:);% Tx + V30
    
    % BMI 5->8
    print_llhd(8,8,:) = mv_logl(5,5,:);% D83 + V30


    mv_logl = print_llhd;
    
      % Print Logl
    figure(fig_num);clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    h_logl = imagesc(mv_logl);
    set(gca,'YDir','normal');
    %set(gca,'XTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});
    %set(gca,'YTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});
    set(gca,'XTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    set(gca,'YTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    set(gca,'FontSize',16);
    set(h_logl,'alphadata',~isnan(mv_logl));
    colorbar;
    %mycmap = get(gcf,'Colormap');
    %set(gcf,'Colormap',flipud(mycmap));
    
    is_nan = isnan(mv_logl(:));
    logl_values = num2str(mv_logl(:),'%6.4g');
    logl_values = strtrim(cellstr(logl_values));%# Remove any space padding
    [x_logl,y_logl] = meshgrid(1:8); %# x/y coordiantes for strings
    h_logl_strings = text(x_logl(~is_nan),y_logl(~is_nan),logl_values(~is_nan),...
            'HorizontalAlignment','center','FontSize',18);
    title('Cox PH Model Log-likelihood','FontSize',18);
    
    
       % Print AIC
%     figure(fig_num);clf reset;
%     fig_num=fig_num+1;
%     set(gcf,'Position',ss_four2three);
%     h_aic = imagesc(mv_aic);
%     set(gca,'YDir','normal');
%     set(gca,'XTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});
%     set(gca,'YTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});    
%     set(gca,'FontSize',16);
%     set(h_aic,'alphadata',~isnan(mv_aic));
%     colorbar;
%     mycmap = get(gcf,'Colormap');
%     set(gcf,'Colormap',flipud(mycmap));
%     
%     is_nan = isnan(mv_aic(:));
%     aic_values = num2str(mv_aic(:),'%6.3g');
%     aic_values = strtrim(cellstr(aic_values));%# Remove any space padding
%     [x_aic,y_aic] = meshgrid(1:8); %# x/y coordiantes for strings
%     h_aic_strings = text(x_aic(~is_nan),y_aic(~is_nan),aic_values(~is_nan),...
%             'HorizontalAlignment','center','FontSize',18);
%     title('Cox PH Model AIC','FontSize',18);
%     
%     
%    
% 	min_info_loss = mv_aic;
%     uni_aics = diag(min_info_loss);
%     
%     prob_min_info_loss = zeros(size(min_info_loss));
%     for i=1:8
%       prob_min_info_loss(i,:) = exp((min_info_loss(i,:)-uni_aics(i))./2);
%     end
%     prob_min_info_loss(prob_min_info_loss>=1)=NaN;
% 
% 
%     % Print AIC Probability Information Loss
%     figure(fig_num);clf reset;
%     set(gcf,'Position',ss_four2three);
%     fig_num=fig_num+1;
%     h_prob_il = imagesc(prob_min_info_loss);
%     set(gca,'YDir','normal');
%     set(gca,'XTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});
%     set(gca,'YTickLabel',{'Tx','Dose/Fx','cm2cw','NumFx','BMI','V_{39}','D_{83}'});  
%     set(gca,'FontSize',16);
%     
%     set(h_prob_il,'alphadata',~isnan(prob_min_info_loss));
%     colorbar;
%     mycmap = get(gcf,'Colormap');
%     set(gcf,'Colormap',flipud(mycmap));
%     
%     is_nan = isnan(prob_min_info_loss(:));
%     il_values = num2str(prob_min_info_loss(:),'%6.2g');
%     il_values = strtrim(cellstr(il_values));%# Remove any space padding
%     [x_il,y_il] = meshgrid(1:8); %# x/y coordiantes for strings
%     h_il_strings = text(x_il(~is_nan),y_il(~is_nan),il_values(~is_nan),...
%             'HorizontalAlignment','center','FontSize',18);
%     title('Relative Prob of Info Loss (compared to univar)','FontSize',18);
%     
    
    % Log-Likelihood Ratio P-values
    
    mv_llr_stats = -2.*mv_logl;
    uni_llr_stats = diag(mv_llr_stats);
    
    llr_pvals1 = zeros(size(mv_llr_stats));
    llr_pvals2 = zeros(size(mv_llr_stats));
    for i=1:8
        cur_llr_stats = uni_llr_stats(i)-mv_llr_stats(i,:);
        llr_pvals1(i,:) = 1-chi2cdf(cur_llr_stats,1); % 1 additional variable
        cur_llr_stats2 = uni_llr_stats(i)-mv_llr_stats(:,i);
        llr_pvals2(:,i) = 1-chi2cdf(cur_llr_stats2,1); % 1 additional variable
    end
    llr_pvals = max(llr_pvals1,llr_pvals2);
    llr_pvals(llr_pvals==1)=NaN;

    
    figure(fig_num);clf reset;
    fig_num=fig_num+1;
    set(gcf,'Position',ss_four2three);
    
    h_llr_pvals = imagesc(llr_pvals);
    set(gca,'YDir','normal');
    set(gca,'XTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    set(gca,'YTickLabel',{'V_{39}','V_{30}','D_{83}','Tx','Dose/Fx','cm2cw','NumFx','BMI'});
    
    set(gca,'FontSize',16);
    set(h_llr_pvals,'alphadata',~isnan(llr_pvals));
    colorbar;
    mycmap = get(gcf,'Colormap');
    set(gcf,'Colormap',flipud(mycmap));
    
    is_nan = isnan(llr_pvals(:));
    llr_pvals(llr_pvals<0.001) = 0.001;
    llr_values = num2str(llr_pvals(:),'%2.2g');
    llr_values = strtrim(cellstr(llr_values));%# Remove any space padding
    [x_il,y_il] = meshgrid(1:8); %# x/y coordiantes for strings
    h_llr_pval_strings = text(x_il(~is_nan),y_il(~is_nan),llr_values(~is_nan),...
            'HorizontalAlignment','center','FontSize',18);
    title('Maximum Log-Likelihood Ratio Statistic p-values','FontSize',18);
    
    %% Tri-variate model
    
     % Tx + cm2cw + bmi
     % Tx + BMI
    [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) cm2cw_data(bmi_idx) bmi_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + cm2cw + BMI');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  cm2cw: ',num2str(tri_beta(2)),...
        '  bmi: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  cm2cw: ',num2str(tri_se(2)),...
        '  bmi: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  cm2cw: ',num2str(tri_p(2)),...
        '  bmi: ',num2str(tri_p(3))]);
    
    
    %% Tx + cm2cw + V_{39}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data cm2cw_data v39_data],...
                compdate,'baseline',0,'censoring',flgcensor);
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + cm2cw + V_{39}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  cm2cw: ',num2str(tri_beta(2)),...
        '  vd: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  cm2cw: ',num2str(tri_se(2)),...
        '  vd: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  cm2cw: ',num2str(tri_p(2)),...
        '  vd: ',num2str(tri_p(3))]);
    
    
    %% Tx + bmi + V_{39}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) bmi_data(bmi_idx) v39_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + bmi + V_{39}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  vd: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  vd: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  vd: ',num2str(tri_p(3))]);
     %% Tx + bmi + V_{30}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) bmi_data(bmi_idx) v30_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + bmi + V_{30}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  v30: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  v30: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  v30: ',num2str(tri_p(3))]);
%% VD + BMI + DosePerFx
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v39_data(bmi_idx) bmi_data(bmi_idx) dosePerFx_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{39} + BMI + DosePerFx');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   V_{39}: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  Dose/Fx: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   V_{39}: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  Dose/Fx: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   V_{39}: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  Dose/Fx: ',num2str(tri_p(3))]);
%% V30 + BMI + DosePerFx
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx) dosePerFx_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{30} + BMI + DosePerFx');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   V_{30}: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  Dose/Fx: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   V_{30}: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  Dose/Fx: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   V_{30}: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  Dose/Fx: ',num2str(tri_p(3))]);
    %% VD + BMI + cm2cw
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v39_data(bmi_idx) bmi_data(bmi_idx) cm2cw_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{39} + BMI + cm2cw');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   V_{39}: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  cm2cw: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   V_{39}: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  cm2cw: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   V_{39}: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  cm2cw: ',num2str(tri_p(3))]);
    %% V30 + BMI + cm2cw
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx) cm2cw_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{30} + BMI + cm2cw');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   V_{30}: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  cm2cw: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   V_{30}: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  cm2cw: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   V_{30}: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  cm2cw: ',num2str(tri_p(3))]);
    
    %% V_{39} + BMI + D_{83}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v39_data(bmi_idx) bmi_data(bmi_idx) d83_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
            
        
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{39} + bmi + D_{83}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   vd: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  dv: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   vd: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  dv: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   vd: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  dv: ',num2str(tri_p(3))]);
    
    %% V_{30} + BMI + D_{83}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx) d83_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
            
        
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('V_{30} + bmi + D_{83}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   v30: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  dv: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   v30: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  dv: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   v30: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  dv: ',num2str(tri_p(3))]);
    
    
     %% Tx + cm2cw + bmi + V_{39}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) cm2cw_data(bmi_idx) bmi_data(bmi_idx) v39_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    quad_logl = cur_logl; 
    quad_aic = -2*cur_logl + 2*4;
    quad_beta = cur_stats.beta;    
    quad_se = cur_stats.se;
    quad_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + cm2cw + bmi + V_{39}');
    disp(['LogL: ',num2str(quad_logl)]);
    disp(['AIC: ',num2str(quad_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(quad_beta(1)),...
        '  cm2cw: ',num2str(quad_beta(2)),...
        '  bmi: ',num2str(quad_beta(3)),...
        '  vd: ',num2str(quad_beta(4))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(quad_se(1)),...
        '  cm2cw: ',num2str(quad_se(2)),...
        '  bmi: ',num2str(quad_se(3)),...
        '  vd: ',num2str(quad_se(4))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(quad_p(1)),...
        '  cm2cw: ',num2str(quad_p(2)),...
        '  bmi: ',num2str(quad_p(3)),...    
        '  vd: ',num2str(quad_p(4))]);
    
    
    
    
     %% V_{30} + BMI + Tx + Dose/Fx
      [~,cur_logl,~,cur_stats]=...
            coxphfit([v30_data(bmi_idx) bmi_data(bmi_idx) tx_data(bmi_idx) dosePerFx_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    quad_logl = cur_logl; 
    quad_aic = -2*cur_logl + 2*4;
    quad_beta = cur_stats.beta;    
    quad_se = cur_stats.se;
    quad_p = cur_stats.p;   
    
%     [quad_p,p_idx] = sort(quad_p);
%     quad_se = quad_se(p_idx);
%     quad_beta = quad_beta(p_idx);
%     quad_aic = quad_aic(p_idx);
%     quad_logl = quad_logl(p_idx);
 
    disp(['']);
    disp('***********');
    disp('V_{30} + BMI + Tx + Dose/Fx');
    disp(['LogL: ',num2str(quad_logl)]);
    disp(['AIC: ',num2str(quad_aic)]);
    disp(['beta: ']);
    disp(['   v30: ',num2str(quad_beta(1)),...
        '  bmi: ',num2str(quad_beta(2)),...
        '  tx: ',num2str(quad_beta(3)),...
        '  Dose/Fx: ',num2str(quad_beta(4))]);
    disp(['se: ']);
    disp(['   v30: ',num2str(quad_se(1)),...
        '  bmi: ',num2str(quad_se(2)),...
        '  tx: ',num2str(quad_se(3)),...
        '  Dose/Fx: ',num2str(quad_se(4))]);
    disp(['p_vals: ']);
    disp(['   v30: ',num2str(quad_p(1)),...
        '  bmi: ',num2str(quad_p(2)),...
        '  tx: ',num2str(quad_p(3)),...    
        '  Dose/Fx: ',num2str(quad_p(4))]);
     
    %% Tx + bmi + V_{99}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) bmi_data(bmi_idx) v99_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + bmi + V_{99}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  v99: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  v99: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  v99: ',num2str(tri_p(3))]);
    
       %% Tx + bmi + V_{165}
      [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data(bmi_idx) bmi_data(bmi_idx) v165_data(bmi_idx)],...
                compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + bmi + V_{165}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  bmi: ',num2str(tri_beta(2)),...
        '  v165: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  bmi: ',num2str(tri_se(2)),...
        '  v165: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  bmi: ',num2str(tri_p(2)),...
        '  v165: ',num2str(tri_p(3))]);
    
    
    %% KPS
       [~,cur_logl,~,cur_stats]=...
            coxphfit([tx_data kps_data v39_data],...
                compdate,'baseline',0,'censoring',flgcensor);
    tri_logl = cur_logl; 
    tri_aic = -2*cur_logl + 2*3;
    tri_beta = cur_stats.beta;    
    tri_se = cur_stats.se;
    tri_p = cur_stats.p;    
 
    disp(['']);
    disp('***********');
    disp('Tx + KPS + V_{39}');
    disp(['LogL: ',num2str(tri_logl)]);
    disp(['AIC: ',num2str(tri_aic)]);
    disp(['beta: ']);
    disp(['   tx: ',num2str(tri_beta(1)),...
        '  kps: ',num2str(tri_beta(2)),...
        '  v39: ',num2str(tri_beta(3))]);
    disp(['se: ']);
    disp(['   tx: ',num2str(tri_se(1)),...
        '  kps: ',num2str(tri_se(2)),...
        '  v39: ',num2str(tri_se(3))]);
    disp(['p_vals: ']);
    disp(['   tx: ',num2str(tri_p(1)),...
        '  kps: ',num2str(tri_p(2)),...
        '  v39: ',num2str(tri_p(3))]);
    
    disp(sprintf('\n'));
    disp(sprintf('\n'));
    disp('     Variable        beta               se            logl          p-val');
    
    disp(uni_print);
    
end


toc;
end