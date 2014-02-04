function ChestWallPainAlpha2BetaLR
tic;

fig_loc = 'Z:/elw/MATLAB/cw_analy/figures/latest/';

%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
%fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

% load data

load(strcat(fp,fn{1}),'CGobj_current');

%CGobj_current.mLymanN = 10.^(-1:0.1:1)';
CGobj_current.mLymanN = 10.^(-2:0.1:0)';

 %% Set dose correction to USCBED, calculate USCBEDs
 for k=1:length([CGobj_current.mGrp])
     CGobj_current.mGrp(k).mBeta2AlphaCorrection = 'BED';
     CGobj_current.mGrp(k).mLymanN = 10.^(-2:0.1:0)';
 end
 
 
a2b = 0.5:0.5:25;
%a2b = 0:0.5:3;

a2b = [Inf a2b];
a2b_llhd = -inf(length(a2b),1);
a2b_loga = -inf(length(a2b),1);
for i=1:length(a2b)
    
    cur_a2b = a2b(i);
    
    CGobj_current.mBeta2Alpha = 1/cur_a2b;
    
    %% Run basic gEUD logistic analysis first
    CGobj_current = CGobj_current.fCalculateEUD();
    % logistic regression analysis
    CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
    
    
    st = [CGobj_current.mLogisticRegressionMat];
    dpf = [st.dev]; % deviations
    st =[st.stats];
    df = [st.dfe]; % degree of freedom
    dpf = dpf./df; % deviations per degree of freedom
    [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
    disp(['best loga for a2b: ',num2str(cur_a2b),' is ',num2str(-log10(CGobj_current.mLymanN(loc)))]);
    loga = -log10(CGobj_current.mLymanN(loc));
   
    
    loglikelyhood = -0.5*dpf;

    
    a2b_llhd(i) =loglikelyhood(loc); % the maximum loglikelyhood

    a2b_loga(i) = loga;
end

fig=figure(1);
set(gcf,'Position',ss_four2three);

a2b_llhd68 = repmat(max(a2b_llhd(2:end))-0.5* 1 /(length(CGobj_current.mGrp)-2),size(a2b));    
    
plot(a2b(2:end),a2b_llhd(2:end),'rs--','LineWidth',2);hold on;
plot(a2b(2:end),a2b_llhd68(2:end),'r--','LineWidth',1);
textbp(['Physical Dose LLHD: ',num2str(a2b_llhd(1),3)],'FontSize',18);
hold off;

end