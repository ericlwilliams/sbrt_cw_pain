function ChestWallPainAlpha2BetaMaxBed

%% Same looping coe as ChestWallPainAlpha2Beta, but runs Cox model on Max BED
% Also, dosen't save to external file but plots after loop, unlike Vd Cox model
% in ChestWallPainAlpha2Beta
tic;

do_print =1;


fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';

%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};

CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

% load data

load(strcat(fp,fn{1}),'CGobj_current');
CGobj = CGobj_current;
CGgrp = CGobj.mGrp;
n_grp = CGobj.mNumInGrp;

%tmp 2011/08/29
a2b = 0:0.1:25;
%a2b = 0:0.5:10;

a2b = [Inf a2b];

% select patients with data
% survival/complication time
f2 = ~cellfun('isempty',{CGgrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CGgrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(n_grp,1);
lastfollowup = inf(n_grp,1);
compdate(f2) = ([CGgrp(f2).mDateComp] - [CGgrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CGgrp(f3).mDateLastFollowup] - [CGgrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [CGgrp.mFlgCensor]';
  
pttotal = ones(CGobj.mNumInGrp,1);
ptcomp = ones(CGobj.mNumInGrp,1); 
ptcomp([CGgrp.mFlgCensor])=0;
            
            
max_bed_cox_logls = zeros(length(a2b),1);
max_bed_lr_logls = zeros(length(a2b),1);

for i=1:length(a2b)
    
    cur_a2b = a2b(i);
    
    max_beds=zeros(n_grp,1);
    %loop of patients, calculate new dosebins and vXs
    for j=1:n_grp
        
        num_fx = CGgrp(j).mFxNum;
        cur_dosebins = CGgrp(j).mDoseBins_org;
        vol = CGgrp(j).mVolCum;
        if cur_a2b<Inf
            new_dosebins = cur_dosebins.* ...
                ((cur_a2b + (cur_dosebins/num_fx))/...
                (cur_a2b + 2));
        else
            new_dosebins = cur_dosebins;
        end
        
        zero_vol = find(vol<=0);
        if isempty(zero_vol) %take last
            cur_max_bed = new_dosebins(end);
        else
            zero_vol_ind = zero_vol(1);
            cur_max_bed = new_dosebins(zero_vol_ind);
        end
        
        max_beds(j)=cur_max_bed;
    end
    %% Cox model
    [~,max_bed_cox_logls(i),~,~] = coxphfit(max_beds,compdate,'baseline',0,'censoring',flgcensor);
    
    %% Logistic Regression
    [~,dev,~]=glmfit(max_beds,[ptcomp pttotal],'binomial','link','logit');
    max_bed_lr_logls(i) = -0.5*dev;
end


cur_fig=figure(1); clf reset; hold on;
set(gcf,'Position',ss_four2three);
plot(a2b(2:end),max_bed_cox_logls(2:end),'bo-');hold on;
textbp(['Phys. Dose Cox LogL: ',num2str(max_bed_cox_logls(1),4),10,...
        '  Low 68% CI: ',num2str(max_bed_cox_logls(1)-0.5,4)],'FontSize',20);
title(['Max BED Cox Model Log-Likelihoods'],'FontSize',22);
set(gca,'FontSize',18);
grid on;
xlabel('LQ model \alpha/\beta','FontSize',22);
ylabel('Cox Model Log-Likelihood','FontSize',22);


  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'cwp_max_bed_a2b_cox_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_max_bed_a2b_cox_llhds.pdf...']);
  end

  
  cur_fig=figure(2); clf reset; hold on;
set(gcf,'Position',ss_four2three);

plot(a2b(2:end),max_bed_lr_logls(2:end),'bo-');hold on;
textbp(['Phys. Dose LogL: ',num2str(max_bed_lr_logls(1),4),10,...
    'Low 68% CI: ',num2str(max_bed_lr_logls(1)-0.5,4)],'FontSize',20);
title(['Max BED Logistic Regression Log-Likelihoods'],'FontSize',22);
set(gca,'FontSize',18);
grid on;
xlabel('LQ model \alpha/\beta','FontSize',22);
ylabel('Logistic Regression Log-Likelihood','FontSize',22);


  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'cwp_max_bed_a2b_lr_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_max_bed_a2b_lr_llhds.pdf...']);
  end


toc;

end