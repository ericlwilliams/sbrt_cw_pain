function CWP_Paper

tic;
fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';
do_print=true;
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

cwp_def = 'MUTTER';

a2b='Inf';
%a2b='2.1';

fig_ctr=1;

%fn = {strcat(cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat')};
fn = {[cwp_def,'_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b',a2b,'.mat']};
%vxdx_cphm_mat_str=strcat(fp,strcat(cwp_def,'_CW_VxDx_CoxPHM.mat'));

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

load(strcat(fp,fn{1}),'CGobj_current');
CGobj = CGobj_current;
clear CGobj_current;

% load basics
dpfx = [CGobj.mGrp.mDosePerFx];
tx_data  = [CGobj.mGrp.mDoseTx]';
fx = [CGobj.mGrp.mFxNum];
flgcensor = [CGobj.mGrp.mFlgCensor]';


[VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
flgCox(flganti)=false; % anti-correlations were not be considered
VDxCox = VDxCox(flgCox);

if a2b=='Inf'
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Histogram of V30
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v30_ind=31;
    if isequal(a2b,'2.1')
        v30_ind=100;
    end
    
    v30_data = [VDxCox(v30_ind).data_exposure];
    
    
    cur_fig=figure(fig_ctr);clf reset;
    set(cur_fig,'Position',ss_four2three);
    fig_ctr=fig_ctr+1;
    hist(v30_data,25);
    disp(['V30Gy = ',num2str(min(v30_data)),' - ',num2str(max(v30_data))]);
    disp(['Median V30Gy = ',num2str(median(v30_data))]);
    set(gca,'FontSize',18);
    xlabel('V_{30Gy} (cc)','FontSize',22);
     
    v30_str = ['Range $V_{30\rm{Gy}} = ',...
        num2str(min(v30_data)),' - ',num2str(round(max(v30_data)*10)/10),'$cc',10,...
        'Median $V_{30\rm{Gy}} = ',num2str(median(v30_data),3),'$cc'];
    text(100,45,v30_str,'FontSize',26,'interpreter','latex');
    
    if do_print
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_loc,'v30'],'-pdf');
        disp(['Saving ',fig_loc,'v30.pdf']);
    end
end %end physical dose plots


% See PlotAlpha2BetaMV.m for bivariate plots/stats


% if a2b=='2.1' % start LQ corrected dose plots
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % V{99} + BMI, CPHM
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%    % Get VDx Cox results     
%     [VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
%     flgCox(flganti)=false; % anti-correlations were not be considered
%     VDxCox = VDxCox(flgCox);
%     logl = [VDxCox.logl]'; 
%     [~,doseloc]=max(logl); % the best fitting of Cox model
%     v99_data = [VDxCox(doseloc).data_exposure];
%     
%     
%     [bmiCox,~,~] = CGobj.fCoxParameter_DVH('BMI');
%     bmi_data = [bmiCox.data_exposure];
%     compdate = [bmiCox.data_hazard];
%     bmi_idx = bmi_data>0;
%     
%     v99_data = v99_data(bmi_idx);
%     bmi_data = bmi_data(bmi_idx);
%     cur_compdate = compdate(bmi_idx);
%     cur_flgcensor = flgcensor(bmi_idx);
%     
%     % V99Gy + BMI CPHM
%     [~,~,~,cur_stats]=...
%             coxphfit([v99_data bmi_data],...
%             cur_compdate,'baseline',0,'censoring',cur_flgcensor); 
%     
%     v99_cph_beta = cur_stats.beta(1);
%     bmi_cph_beta = cur_stats.beta(2);
%         
%     disp(['V99_{Gy_{2.1}} + BMI CPHM Results']);
%     disp(['Beta (V99, BMI): (',...
%         num2str(v99_cph_beta),',',...
%         num2str(bmi_cph_beta),')']);
%     disp(['P (V99, BMI): (',...
%         num2str(cur_stats.p(1)),',',...
%         num2str(cur_stats.p(2)),')']);
%      
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % V{99} + BMI, Logrank
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     vd_bmi = v99_cph_beta.*v99_data +...
%         bmi_cph_beta.*bmi_data;
%     
%     
%     vd_bmi_split = median(vd_bmi);
%     
%     below_median= (vd_bmi<vd_bmi_split);
%     
%   
%     % initialize a survivalanalysis obj
%     sa=classKaplanMeierCurve(); 
%     survivedate={cur_compdate(below_median); cur_compdate(~below_median)}; % survive time of each group
%     fcensor={cur_flgcensor(below_median); cur_flgcensor(~below_median)}; % censor flag for each group
%     sa.mSurvivalTime=survivedate;
%     sa.mFlgCensor=fcensor;
%     
%     % compute survival curves and compare them
%     sa=sa.fCalculateSurvivalCurve();
%     sa=sa.fCombineSurvivalTime();
%     sa=sa.fCompareSurvivalByLogrank();
%    
%     cox_beta=coxphfit(~below_median,cur_compdate,'baseline',0,'censoring',cur_flgcensor);
%     cox_hr = exp(cox_beta);
%     
%     vd_bmi_lr_pval = sa.mpValue;
%     disp(['HR: ',num2str(cox_hr)]);
%     disp(['Logrank p: ',num2str(vd_bmi_lr_pval)]);
%     
%     %plot split
%     figure(fig_ctr);  clf reset; hold on; % grid on;
%     fig_ctr=fig_ctr+1;
%     h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
%     plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
%         1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
%     h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
%     plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
%         1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
%     %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
%     text(38,0.25,str_pval1,'FontSize',12);
%     lgnd=legend(h_km,...
%         strcat('V$_{39}\leq',num2str(vol,4),'$'),...
%         strcat('V$_{39}\geq',num2str(vol,4),'$'),'Location','Best');
%     set(lgnd,'FontSize',14);
%     h=legend;
%     set(h,'interpreter','latex');
%     
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel(['Months'],'fontsize',14);
%     ylabel(['Probability of CW Pain'],'fontsize',14);
%     title('V_{39}, Median split','fontsize',14);
%     
% end

end