function ChestWallUscDtProfile
tic;
% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

do_print = 1;

fig_loc = 'Z:\elw\MATLAB\cw_analy\slides\figures\latest\';
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
%fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'};
CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
load(strcat(fp,fn{1}),'CGobj_current');
%CGgrp = CGobj_current.mGrp;

CGobj_current.mLymanN = 10.^(-1:0.1:1)';

%LymanN = log10(CGmsk.mLymanN);
%CGobj_current.LymanN = log10(CGmsk.mLymanN);

% %% Set dose correction to USCBED, calculate USCBEDs
 for k=1:length([CGobj_current.mGrp])
     CGobj_current.mGrp(k).mBeta2AlphaCorrection = 'USCBED';
     CGobj_current.mGrp(k).mLymanN = 10.^(-1:0.1:1)';
 end
 
 
 %% testing
 
 %mUscRangeDt = 0.1:1:30;
 %mUscRangeAlphaD0 = 0.01:0.1:.5;
 
 
 mUscRangeDt = 0.1:0.5:30;
 mUscRangeAlphaD0 = 0.01:0.05:.5;
 
 
 
 mUscAlpha2Beta = 3; 

mUscBestLog10a = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscLogLikelihoods = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscLogLikelihoods68 = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscLogLikelihoods95 = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscAICs = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscDq = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscPvals = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscFracLQ = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));
mUscFracFullLQ = inf(length(mUscRangeDt),length(mUscRangeAlphaD0));

 for i=1:length(mUscRangeDt)
        
     cur_dt = mUscRangeDt(i);
        
     for j=1:length(mUscRangeAlphaD0)
        
        cur_ad0 = mUscRangeAlphaD0(j);

        cur_dq = 0.5*cur_dt*(1-cur_ad0);
        
        [CGobj_current,cur_frac_lq,cur_frac_full_lq] = CGobj_current.fUSCCorrection(cur_dq,cur_ad0,mUscAlpha2Beta);
        
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
        %
        
        [mxllhd,~] = max(loglikelyhood); % the maximum loglikelyhood
        
        loglikelyhood68 = mxllhd-0.5* 1 /df(loc);
        loglikelyhood95 = mxllhd-0.5* (1.96*2) /df(loc);
        
        
        aic = -2*loglikelyhood.*df + 2*3; % 3 parameters in usc
        [minaic,~] = min(aic);
        %disp(['Max LLHD: ',num2str(mxllhd),...
        %    ' for d0: ',num2str(cur_d0),...
        %    ' dq: ',num2str(cur_dq),...
        %    ' dt: ',num2str(cur_dt)]);
        
        mUscBestLog10a(i,j) = loga;
        mUscLogLikelihoods(i,j) = mxllhd;
        mUscLogLikelihoods68(i,j) = loglikelyhood68;
        mUscLogLikelihoods95(i,j) = loglikelyhood95;
        mUscAICs(i,j) = minaic;
        mUscDq(i,j) = cur_dq;
        mUscPvals(i,j) = min_p;
        mUscFracLQ(i,j) = cur_frac_lq;
        mUscFracFullLQ(i,j) = cur_frac_full_lq;
         
     end
 end
    fig_ctr=0;
   fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
        
    llhds = max(mUscLogLikelihoods,[],2)';
    llhd68 = max(max(mUscLogLikelihoods68,[],2));
    llhd95 = max(max(mUscLogLikelihoods95,[],2));    
    plot(mUscRangeDt,llhds,'--','LineWidth',2);hold on;
    h_llhd68=plot(mUscRangeDt,repmat(llhd68,1,length(mUscRangeDt)),'g--','LineWidth',1);
    h_llhd95=plot(mUscRangeDt,repmat(llhd95,1,length(mUscRangeDt)),'r--','LineWidth',1);
    hold off;
    
    %interp1
    %plot(mUscRangeAlphaD0,mUscLogLikelihoods,[],2)','--','LineWidth',2);hold on;
    
    llhd_lgnd = legend([h_llhd68 h_llhd95],'68% CI','95% CI','location','best');
    set(llhd_lgnd,'FontSize',18);
    set(gca,'FontSize',18);
    xlabel('D_{T} [Gy]','FontSize',22);
    ylabel('Max Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_prof_mxllhd_vs_dt'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_prof_mxllhd_vs_dt.pdf...']);
  end
 
end
 