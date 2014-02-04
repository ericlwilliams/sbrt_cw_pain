function ChestWallCoxFitDisplay
tic;

% prepare
%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\Desktop\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';


fig_loc = 'Z:\elw\MATLAB\cw_analy\figures\latest\';
%fn = {'MUTTER_MASTER_VDx_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
%fn = {'MUTTER_MASTER_VDx_ChestWall_Cox_DiVj_DVHs_fx-1_a2b3.mat'};
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
CGobj = cell(length(fn),1);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end


for m = 1:length(fn)
    OCobj=CGobj{m};
    [allCox,flgCox,flganti] = OCobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
    flgCox(flganti)=false; % anti-correlations were not be considered
    logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
    [~,doseloc]=max(logl); % the best fitting of Cox model
    allCox = allCox(doseloc);
    disp(['Best Cox Model is at dose :',num2str(OCobj.mBinsDose(doseloc))]);
disp(allCox);
    
    % Get Vx at best dose
    % volumes of patients at best Vx
    Vx=zeros(OCobj.mNumInGrp,1);
    d = OCobj.mBinsDose(doseloc);
    for k=1:OCobj.mNumInGrp
        Vx(k) = OCobj.mGrp(k).fVolAtDose( d );
    end
    
    hzrd_fn = allCox.h;
    
    %adjust duplicate values of t
    f = find(diff(hzrd_fn(:,1))==0); % find duplicate time values of h(t)
    while ~isempty(f)
        hzrd_fn(f,1) = hzrd_fn(f,1)-eps*10; % adjust it a bit to avoid ambiguius
        f = find(diff(hzrd_fn(:,1))==0); % find duplicate time values of h(t)
    end
            
    
    cox_t = hzrd_fn(:,1);
    cox_h = hzrd_fn(:,2);
    
     % survival curve from Cox model at volume vol
    km_t = OCobj.mKaplanMeierCompOverall.mSurvivalTimeSorted{1};
    cens_idx=OCobj.mKaplanMeierCompOverall.mCensorStatistics{1}(:,1);
    km_t(cens_idx) = []; %remove censored times, left with comp times
    
    % km_t 47x1 because from complication frequency table (means 47 unique
    % comp times? see mSurvialTimeTable in KaplanMeier
    cox_comp_h = interp1(cox_t,cox_h,km_t,'linear','extrap');
    
    
    S0 = exp(-cox_comp_h);
    cox_s = S0.^(exp(allCox.beta*mean(Vx)));
    
   % expbetax = exp(CoxPar.beta*vol);
    %CoxComplicationCurve = exp( -h * expbetax );
    
    %% KM curve
    
    % complication incidence curve
disp('Complication incidence curves');
m = 1;
figure(m); clf reset; hold on; % grid on;
sa = OCobj.mKaplanMeierCompOverall;
stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');

% stairs(km_t,1-sa.mSurvivalCurve{1});
% plot(km_t(sa.mCensorStatistics{1}(:,1)),...
%     1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');

plot(km_t,1-cox_s,'r--');
%         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
set(gca,'xminortick','on','yminortick','on');
xlabel(['Months']);
ylabel(['Probability of grade >= 2 Chestwall Pain']);


    
end

