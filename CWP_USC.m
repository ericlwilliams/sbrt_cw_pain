function CWP_USC
tic;
% prepare
fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
%fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

do_print = 1;

fig_loc = 'Z:\elw\MATLAB\cw_analy\slides\figures\latest\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};

%fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'};
CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
load(strcat(fp,fn{1}),'CGobj_current');
%CGgrp = CGobj_current.mGrp;

do_lq=false;

CGobj_current.mLymanN = 10.^(-2:0.1:0)';

%LymanN = log10(CGmsk.mLymanN);
%CGobj_current.LymanN = log10(CGmsk.mLymanN);

% %% Set dose correction to USCBED, calculate USCBEDs
 for k=1:length([CGobj_current.mGrp])
     CGobj_current.mGrp(k).mBeta2AlphaCorrection = 'BED';
     CGobj_current.mGrp(k).mLymanN = 10.^(-2:0.1:0)';
 end
%% Physical gEUD
%% Run basic gEUD logistic analysis first

CGobj_phys = CGobj_current;
 CGobj_phys = CGobj_phys.fCalculateEUD();
        % logistic regression analysis
 CGobj_phys = CGobj_phys.fLogisticRegressionExact_EUD();

 %% printing lq
 CGobj_phys.mLymanN = log10(CGobj_phys.mLymanN);
 
 fig_ctr=1;
 cur_fig=figure(fig_ctr);  clf reset;% grid on;
 set(gcf,'Position',ss_four2three);

 [loga,pval]=CGobj_phys.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','r--',2); 
    CGobj_phys.fComplicationObservedFig_EUD(loga,4,'k*',2);
    ylim([0 0.5]);
xlabel(['gEUD_{PHYS} log_{10}(a) =',num2str(loga,2)],'FontSize',22);
ylim([0 0.5]);
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'phys_geud_resp'],'-pdf');
   disp(['Saving ',fig_loc,'phys_geud_resp.pdf...']);
end


%% LQ gEUD
if do_lq
 CGobj_current.mBeta2Alpha = 1/3;
end
 %% Run basic gEUD logistic analysis first
 CGobj_lq = CGobj_current.fCalculateEUD();
        % logistic regression analysis
 CGobj_lq = CGobj_lq.fLogisticRegressionExact_EUD();

 %% printing lq
 CGobj_lq.mLymanN = log10(CGobj_lq.mLymanN);
 
 fig_ctr=1;
 cur_fig=figure(fig_ctr);  clf reset;% grid on;
 set(gcf,'Position',ss_four2three);
[~, ~,lq_geud_llhd,lq_geud_loga] = CGobj_lq.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','rs--',2);
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_llhds.pdf...']);
 end

 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
[loga,pval]=CGobj_lq.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','r--',2); 
CGobj_lq.fComplicationObservedFig_EUD(loga,4,'k*',2);
ylim([0 0.5]);
xlabel(['gEUD_{LQ} log_{10}(a) =',num2str(loga,2)],'FontSize',22);
ylim([0 0.5]);
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_resp'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_resp.pdf...']);
end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
    set(gcf,'Position',ss_four2three);
 CGobj_lq.fLogisticRegressionPvalueExactFig_a_EUD('rs--',2);
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_pvals'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_pvals.pdf...']);
end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 [~,lq_aic_llhd,lq_geud_aic,lq_aic_log10a] =  CGobj_lq.fLogisticRegressionAicFig_a_EUD('loga',1,'rs--',2); % 1 parameter for LQ model
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_aics'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_aics.pdf...']);
 end
 
 %% Max BED
 pttotal = ones(CGobj_lq.mNumInGrp,1);
 ptcomp = ones(CGobj_lq.mNumInGrp,1);
 CGgrp = [CGobj_lq.mGrp];
 
 ptcomp([CGgrp.mFlgCensor])=0;
 
 max_beds = -inf(length(CGgrp),1);
 for j=1:length(CGgrp)
     cur_dosebins = CGgrp(j).mDoseBins_LQ;
     vols = CGgrp(j).mVolCum;
     
     zero_vol = find(vols<=0);
     if isempty(zero_vol) %take last
         cur_max_bed = cur_dosebins(end);
     else
         zero_vol_ind = zero_vol(1);
         cur_max_bed = cur_dosebins(zero_vol_ind);
     end
     max_beds(j) = cur_max_bed;
 end
 
 %% Logistic Regression
 [b,dev,st]=glmfit(max_beds,[ptcomp pttotal],'binomial','link','logit');
 
 max_bed_lr_logl = -0.5*dev/(st.dfe);
 disp(['Max BED LLHD: ',num2str(max_bed_lr_logl,3)]);
 
 lq_dm_llhd = max_bed_lr_logl;
 lq_dm_aic = -2*lq_dm_llhd*st.dfe+2;
 doses = (0:max(max_beds))'; % doses (gEUD) whose RP probability will be computed
 
 [rpb,rplo,rphi] = glmval(b, doses,'logit',st); % the responding function values at doses
 
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 plot(doses,rpb,'r--','LineWidth',2);   hold on; % responding function
 plot(doses,rpb-rplo,'r--','LineWidth',1); % low CI curve
 plot(doses,rpb+rphi,'r--','LineWidth',1); % high CI curve
 
 flg=[CGobj_lq.mGrp.mFlgCensor]; % censor flags of patients
 [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,max_beds,4);
 prob = numcomp./numtotal;
 % plot
 errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',2);
 errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
 hold off;
 ylim([0 0.5]);
 set(gca,'xminortick','on','yminortick','on');
 set(gca,'box','on');
 set(gca,'FontSize',18);
 xlabel('Max BED Dose [Gy]','FontSize',20); ylabel('CWP probability','FontSize',20);
   
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_dmax_resp'],'-pdf');
   disp(['Saving ',fig_loc,'lq_dmax_resp.pdf...']);
 end
 
 
 
 %% Max phys
  
 max_phys = -inf(length(CGgrp),1);
 for j=1:length(CGgrp)
     cur_dosebins = CGgrp(j).mDoseBins_org;
     vols = CGgrp(j).mVolCum;
     
     zero_vol = find(vols<=0);
     if isempty(zero_vol) %take last
         cur_max_phys = cur_dosebins(end);
     else
         zero_vol_ind = zero_vol(1);
         cur_max_phys = cur_dosebins(zero_vol_ind);
     end
     max_phys(j) = cur_max_phys;
 end
 
 %% Logistic Regression
 [b,dev,st]=glmfit(max_phys,[ptcomp pttotal],'binomial','link','logit');
 
 max_phys_lr_logl = -0.5*dev/(st.dfe);
 disp(['Max PHYS LLHD: ',num2str(max_phys_lr_logl,3)]);
 
 doses = (0:max(max_phys))'; % doses (gEUD) whose RP probability will be computed
 
 [rpb,rplo,rphi] = glmval(b, doses,'logit',st); % the responding function values at doses
 
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 plot(doses,rpb,'r--','LineWidth',2);   hold on; % responding function
 plot(doses,rpb-rplo,'r--','LineWidth',1); % low CI curve
 plot(doses,rpb+rphi,'r--','LineWidth',1); % high CI curve
 
 [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,max_phys,4);
 prob = numcomp./numtotal;
 % plot
 errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',2);
 errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
 hold off;
 ylim([0 0.5]);
 set(gca,'xminortick','on','yminortick','on');
 set(gca,'box','on');
 set(gca,'FontSize',18);
 xlabel('Max PHYS','FontSize',20); ylabel('CWP probability','FontSize',20);
   
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'phys_dmax_resp'],'-pdf');
   disp(['Saving ',fig_loc,'phys_dmax_resp.pdf...']);
 end
 %set alpha/beta AND calculate all USC BEDs
 %CGobj_current.mBeta2Alpha(mUscAlpha);
 
 %mUscAlpha = 0.052; % using 0rib cage, see http://www.rooj.com/Normal%20Tissue%20Comp.htm
 mUscAlpha2Beta = 3;
 
 %testing
 %mUscRangeAlphaD0 = 0.01:0.01:0.1;
%mUscRangeDq = 0.2:0.2:1;

%best CWP range
mUscRangeAlphaD0 = 0.21:0.01:0.25;
mUscRangeDq = 6.2:0.2:6.8;


%best Bpx range
%mUscRangeAlphaD0 = 0.01:0.01:0.05;
%mUscRangeDq = 4.5:0.2:6.0;

% normal
%mUscRangeAlphaD0 = 0.01:0.01:0.5;
%mUscRangeDq = 0.2:0.2:7.6;


%% Loop over D0 Dq values, calculate DT, compute USC BED, fit outcome data           
mUscBestLog10a = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscLogLikelihoods = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscAICs = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscDt = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscPvals = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscFracLQ = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscFracFullLQ = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));


mUscDmaxLogLikelihoods = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscDmaxAICs = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscDmaxPvals = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));

for i=1:length(mUscRangeAlphaD0)
    cur_alphad0 = mUscRangeAlphaD0(i);
    for j=1:length(mUscRangeDq)
        cur_dq = mUscRangeDq(j);
        
        % calc USC BEDs
        [CGobj_current,cur_frac_lq,cur_frac_full_lq] = CGobj_current.fUSCCorrection(cur_dq,cur_alphad0,mUscAlpha2Beta);
        
        %% Dmax
        CGobj_dmax = CGobj_current; % don't do gEUD calc, see below
        
        %shortcut
        usc_dmax = [max([CGobj_dmax.mGrp.mDoseBins_LQ])];
 
    %% Logistic Regression
        [~,dev,st]=glmfit(usc_dmax,[ptcomp pttotal],'binomial','link','logit');
                
        mUscDmaxLogLikelihoods(i,j) = -0.5*dev/(st.dfe);
        mUscDmaxAICs(i,j) = -2*mUscDmaxLogLikelihoods(i,j).*(st.dfe) + 2*3; % 3 parameters in usc
        mUscDmaxPvals(i,j) = st.p(2);       
        
        %% gEUD analysis
        % calc gEUD
        CGobj_current = CGobj_current.fCalculateEUD();
        % logistic regression analysis
        CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
        
        %% Find best log10a and LLHDs
        st = [CGobj_current.mLogisticRegressionMat];
        dpf = [st.dev]; % deviations
        st =[st.stats];
        pvalue = [st.p];
        pvalue = pvalue(2,:); % the p-value corresponding to gEUD
        [min_p,~] = min(pvalue);
        
        df = [st.dfe]; % degree of freedom
        dpf = dpf./df; % deviations per degree of freedom
        [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
        %disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-CGobj_current.mLymanN(loc))]);
        loga = -log10(CGobj_current.mLymanN(loc));
        
        %disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
        
        loglikelyhood = -0.5*dpf;
        [mxllhd,~] = max(loglikelyhood); % the maximum loglikelyhood
        
        aic = -2*loglikelyhood.*df + 2*3; % 3 parameters in usc
        [minaic,~] = min(aic);
        %disp(['Max LLHD: ',num2str(mxllhd),...
        %    ' for d0: ',num2str(cur_d0),...
        %    ' dq: ',num2str(cur_dq),...
        %    ' dt: ',num2str(cur_dt)]);
        
        mUscBestLog10a(i,j) = loga;
        mUscLogLikelihoods(i,j) = mxllhd;
        mUscAICs(i,j) = minaic;
        mUscDt(i,j) = (2*cur_dq)/(1-(cur_alphad0));
        mUscPvals(i,j) = min_p;
        mUscFracLQ(i,j) = cur_frac_lq;
        mUscFracFullLQ(i,j) = cur_frac_full_lq;
        
    
    end
    
    % plotyy as function of mUscDt?  will have multiple entries for each
    % mUscDt?
end
[best_llhd, mx_ind] = max(mUscLogLikelihoods(:));
[alphad0_ind, dq_ind] = ind2sub(size(mUscLogLikelihoods),mx_ind);

best_dt = mUscDt(alphad0_ind,dq_ind);
best_pval = mUscPvals(alphad0_ind,dq_ind);
best_log10a = mUscBestLog10a(alphad0_ind,dq_ind);
best_aic = mUscAICs(alphad0_ind,dq_ind);
best_alphad0 = mUscRangeAlphaD0(alphad0_ind);
best_dq = mUscRangeDq(dq_ind);
disp([]);
disp([]);
disp(['== Best USC model ==',10,...
    'LLHD: ',num2str(best_llhd),10,...
     'AIC: ',num2str(best_aic),10,...
    'log10a: ',num2str(best_log10a),10,...
    'pval: ',num2str(best_pval),10,...
    'alpha*D0: ',num2str(best_alphad0),10,...
    'Dq: ',num2str(best_dq),10,...
    'Dt: ',num2str(best_dt)]);
disp([]);
disp([]);


%% USC Dmax
[best_dm_llhd, mx_dm_ind] = max(mUscDmaxLogLikelihoods(:));
[alphad0_dm_ind, dq_dm_ind] = ind2sub(size(mUscDmaxLogLikelihoods),mx_dm_ind);

best_dm_dt = mUscDt(alphad0_dm_ind,dq_dm_ind);
best_dm_pval = mUscPvals(alphad0_dm_ind,dq_dm_ind);
best_dm_alphad0 = mUscRangeAlphaD0(alphad0_dm_ind);
best_dm_dq = mUscRangeDq(dq_dm_ind);
best_dm_aic = mUscDmaxAICs(alphad0_dm_ind,dq_dm_ind);
disp([]);
disp([]);
disp(['== Best USC Dmax model ==',10,...
    'LLHD: ',num2str(best_dm_llhd),10,...
    'AIC: ',num2str(best_dm_aic),10,...
        'pval: ',num2str(best_dm_pval),10,...
    'alpha*D0: ',num2str(best_dm_alphad0),10,...
    'Dq: ',num2str(best_dm_dq),10,...
    'Dt: ',num2str(best_dm_dt)]);
disp([]);
disp([]);




[best_aic, min_aic_ind] = min(mUscAICs(:));
[alphad0_aic_ind, dq_aic_ind] = ind2sub(size(mUscAICs),min_aic_ind);
best_aic_llhd = mUscLogLikelihoods(alphad0_aic_ind,dq_aic_ind);
best_aic_dt = mUscDt(alphad0_aic_ind,dq_aic_ind);
best_aic_pval = mUscPvals(alphad0_aic_ind,dq_aic_ind);
best_aic_log10a = mUscBestLog10a(alphad0_aic_ind,dq_aic_ind);
best_aic_alphad0 = mUscRangeAlphaD0(alphad0_aic_ind);
best_aic_dq = mUscRangeDq(dq_aic_ind);

    
disp(['== AICs == ',10,...
    'USC min AIC: ',num2str(best_aic),10,...
     '- LLHD: ',num2str(best_aic_llhd),10,...
    '- log10a: ',num2str(best_aic_log10a),10,...
    '- pval: ',num2str(best_aic_pval),10,...
    '- alpha*D0: ',num2str(best_aic_alphad0),10,...
    '- Dq: ',num2str(best_aic_dq),10,...
    '- Dt: ',num2str(best_aic_dt)]);
disp([]);
disp(['LQ min AIC: ',num2str(lq_geud_aic),10,...
    '- log10a: ',num2str(lq_aic_log10a),10,...
    '- LLHD: ',num2str(lq_aic_llhd)]);
        disp([]);
            disp([]);

fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscDmaxLogLikelihoods);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Log-likelihood/df','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
title('USC Dmax');
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_dmax_geud_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'usc_dmax_geud_llhds.pdf...']);
 end

 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscLogLikelihoods);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Log-likelihood/df','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_llhds.pdf...']);
 end

 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscFracLQ);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Fraction LQ Bins','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_fraclq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_fraclq.pdf...']);
 end
 
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscFracFullLQ);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Fraction Full CW LQ ','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_fracfulllq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_fracfulllq.pdf...']);
 end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscAICs);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'AIC','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_aics'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_aics.pdf...']);
 end
  
 
fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscBestLog10a);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Best log10(a)','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_log10a'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_log10a.pdf...']);
  end
toc;

fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscDt);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'DT','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_dt'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_dt.pdf...']);
 end
  
 
fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscPvals);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'P-value','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_pvals.pdf...']);
 end
  
 %% get best USC response
 % calc USC BEDs
 [CGobj_current,best_frac_lq,best_frac_full_lq] = CGobj_current.fUSCCorrection(best_dq,best_alphad0,mUscAlpha2Beta);
        
    % calc gEUD
 CGobj_current = CGobj_current.fCalculateEUD();
        % logistic regression analysis
  CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
 
  fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
 CGobj_current.mLymanN = log10(CGobj_current.mLymanN);
    [~,pval]=CGobj_current.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','r--',2); 
    CGobj_current.fComplicationObservedFig_EUD(best_log10a,4,'k*',2);
    
    ylim([0 0.5]);
    xlabel(['gEUD_{USC} log_{10}(a) =',num2str(best_log10a,1)],'FontSize',22);

    disp([]);
disp(['Best USC response',10,...
    'p = ',num2str(pval),10,...
    'log10a = ',num2str(best_log10a),10,...
    'Frac LQ Bins= ',num2str(best_frac_lq),10,...
    'Frac Full CW LQ = ',num2str(best_frac_full_lq)]);
    disp([]);
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'usc_geud_resp'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_resp.pdf...']);
 end
 
%% get best USC Dmax response
 % calc USC BEDs
 [CGobj_current,~,~] = CGobj_current.fUSCCorrection(best_dm_dq,best_dm_alphad0,mUscAlpha2Beta);
 
 
  usc_dmax = -inf(length([CGobj_current.mGrp]),1);
    for j=1:length(CGobj_current.mGrp)
            cur_dosebins = CGobj_current.mGrp(j).mDoseBins_LQ;
            vols = CGobj_current.mGrp(j).mVolCum;
     
            zero_vol = find(vols<=0);
           if isempty(zero_vol) %take last
             cur_usc_dmax = cur_dosebins(end);
           else
            zero_vol_ind = zero_vol(1);
            cur_usc_dmax = cur_dosebins(zero_vol_ind);
           end
        usc_dmax(j) = cur_usc_dmax;
        end
 
 %% Logistic Regression
 [b,~,st]=glmfit(usc_dmax,[ptcomp pttotal],'binomial','link','logit');
 
 doses = (0:max(usc_dmax))'; % doses (gEUD) whose RP probability will be computed
 
 [rpb,rplo,rphi] = glmval(b, doses,'logit',st); % the responding function values at doses
 
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 plot(doses,rpb,'r--','LineWidth',2);   hold on; % responding function
 plot(doses,rpb-rplo,'r--','LineWidth',1); % low CI curve
 plot(doses,rpb+rphi,'r--','LineWidth',1); % high CI curve
 
 flg=[CGobj_lq.mGrp.mFlgCensor]; % censor flags of patients
 [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,usc_dmax,4);
 prob = numcomp./numtotal;
 % plot
 errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',2);
 errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
 hold off;
 
 set(gca,'xminortick','on','yminortick','on');
 set(gca,'box','on');
 set(gca,'FontSize',18);
 xlabel('USC Dmax [Gy]','FontSize',20); ylabel('CWP probability','FontSize',20);
   
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'usc_dmax_resp'],'-pdf');
   disp(['Saving ',fig_loc,'usc_dmax_resp.pdf...']);
 end
 
 
 
 %% profile dt
 
 dt_range = [min(mUscDt):mean(mean(diff(mUscDt))):max(mUscDt(:))]';
 dt_llhds = inf(length(dt_range),1);
 dt_pvals = inf(length(dt_range),1);
 dt_fraclq = inf(length(dt_range),1);
 dt_fracfulllq = inf(length(dt_range),1); 
 for i=1:length(dt_range)
    cur_dt = dt_range(i);
    
    [~,dt_ind] = min(abs(mUscDt(:) - cur_dt));
     
     [alphad0_dt_ind, dq_dt_ind] = ind2sub(size(mUscDt),dt_ind);
      %dt_alphad0 = mUscRangeAlphaD0(alphad0_dt_ind);     
      %dt_dq = mUscRangeDq(dq_dt_ind);     
      dt_llhds(i) = mUscLogLikelihoods(alphad0_dt_ind,dq_dt_ind);
      dt_pvals(i) = mUscPvals(alphad0_dt_ind,dq_dt_ind);
      dt_fraclq(i) = mUscFracLQ(alphad0_dt_ind,dq_dt_ind);
      dt_fracfulllq(i) = mUscFracFullLQ(alphad0_dt_ind,dq_dt_ind);
 end
 
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
 
    plot(dt_range,dt_llhds,'--','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('D_{T} [Gy]','FontSize',22);
    ylabel('Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_dt'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_dt.pdf...']);
  end
  
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
    
    plot(mUscRangeDq,max(mUscLogLikelihoods,[],1),'--','LineWidth',2);
    
    set(gca,'FontSize',18);
    xlabel('D_{q} [Gy]','FontSize',22);
    ylabel('Max Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_dq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_dq.pdf...']);
  end
  
    
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
    
    plot(mUscRangeAlphaD0,max(mUscLogLikelihoods,[],2)','--','LineWidth',2);
    
    set(gca,'FontSize',18);
    xlabel('\alpha\cdot D_{0} [Gy]','FontSize',22);
    ylabel('Max Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_alphad0'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_alphad0.pdf...']);
  end
  
  fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
 
    h_full_frac_lq = plot(dt_range,dt_fracfulllq,'b--','LineWidth',2);hold on;
    h_full_bin_lq = plot(dt_range,dt_fraclq,'g--','LineWidth',2);
    hold off;
    ylim([0 1]);
   xlabel('D_{T} [Gy]','FontSize',22);

   ylabel('Fraction LQ','FontSize',22);
   lgnd = legend([h_full_frac_lq, h_full_bin_lq],'Fraction Full CW LQ','Fraction Bins LQ','location','best');
   set(lgnd,'FontSize',18);
   set(gca,'FontSize',18);
   
  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_lq_fracs_vs_dt'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_lq_fracs_vs_dt.pdf...']);
  end
 
 disp(['*** Beamer Vars ***']);
 if do_lq
    disp(['Using LQ doses']);
 else
     disp(['Using PHYS doses']);
 end

 %% Copy/paste into model_params.tex
 disp([]);
 disp(['\newcommand{\UscGeudLlhd}{',num2str(best_llhd,'%4.3f'),'}']);
disp(['\newcommand{\UscGeudLoga}{',num2str(best_log10a,'%2.1f'),'}']);
disp(['\newcommand{\UscGeudAIC}{',num2str(best_aic,'%5.2f'),'}']);
disp(['\newcommand{\UscGeudAlphaDzero}{',num2str(best_alphad0,'%3.2f'),'}']);
disp(['\newcommand{\UscGeudDq}{',num2str(best_dq,'%3.1f'),'}']);
disp(['\newcommand{\UscGeudDt}{',num2str(best_dt,'%3.1f'),'}']);
 disp([]);
  
 disp(['\newcommand{\UscDmaxLlhd}{',num2str(best_dm_llhd,'%4.3f'),'}']);
disp(['\newcommand{\UscDmaxAIC}{',num2str(best_dm_aic,'%5.2f'),'}']);
disp(['\newcommand{\UscDmaxAlphaDzero}{',num2str(best_dm_alphad0,'%3.2f'),'}']);
disp(['\newcommand{\UscDmaxDq}{',num2str(best_dm_dq,'%3.1f'),'}']);
disp(['\newcommand{\UscDmaxDt}{',num2str(best_dm_dt,'%3.1f'),'}']);
 

if ~do_lq
    disp('= Phys =');
    disp(['\newcommand{\PhysGeudLlhd}{',num2str(lq_geud_llhd,'%4.3f'),'}']);
    disp(['\newcommand{\PhysGeudLoga}{',num2str(lq_geud_loga,'%2.1f'),'}']);
    disp(['\newcommand{\PhysGeudAIC}{',num2str(lq_geud_aic,'%5.2f'),'}']);
    
    disp(['\newcommand{\PhysDmaxLlhd}{',num2str(lq_dm_llhd,'%4.3f'),'}']);
     disp(['\newcommand{\PhysDmaxAIC}{',num2str(lq_dm_aic,'%5.2f'),'}']);
     disp([]);
     disp(['Phys gEUD and Dmax only, need to run LQ']);
else
    disp('= LQ =');
    disp(['\newcommand{\LqGeudLlhd}{',num2str(lq_geud_llhd,'%4.3f'),'}']);
    disp(['\newcommand{\LqGeudLog10a}{',num2str(lq_geud_loga,'%2.1f'),'}']);
    disp(['\newcommand{\LqGeudAIC}{',num2str(lq_geud_aic,'%5.2f'),'}']);
    
    disp(['\newcommand{\LqDmaxLlhd}{',num2str(lq_dm_llhd,'%4.3f'),'}']);
    disp(['\newcommand{\LqDmaxAIC}{',num2str(lq_dm_aic,'%5.2f'),'}']);
    
    disp(['LQ gEUD and Dmax only, need to run Phys']);
end
toc;

end