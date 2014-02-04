function CWP_KMphys
tic;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

fig_loc = 'Z:/elw/MATLAB/cw_analy/slides/figures/latest/';

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

%fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'];

fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'];
%fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx4_a2b2.1.mat'];



load(strcat(fp,fn2),'CGobj_current');
CGobj = {CGobj_current};
clear CGobj_current;
CGobj = CGobj{1};

% Find LR, KM, and HR for given Vd and fraction size

% select patients with data
f = CGobj.fPatientsWithComplicationData();
CGobj = CGobj.fRemovePatient(~f);

% survival/complication time
f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date

compdate = inf(CGobj.mNumInGrp,1);
lastfollowup = inf(CGobj.mNumInGrp,1);
compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [CGobj.mGrp.mFlgCensor]';


numFx = [3 4 5];
%doses = [32 36 40]; %from V99 a2b=2.1
%numFx =[4];
doses = [99 99 99]; %from V99 a2b=2.1
%doses = [99]; %from V99 a2b=2.1
phys_doses = [32 36 40];

grp =[CGobj.mGrp];
fxs = [grp.mFxNum];

for i=1:length(numFx)
    cur_fxs = fxs==numFx(i);
    cur_grp = grp(cur_fxs);
    
    cur_compdate = compdate(cur_fxs);
     cur_flgcensor = flgcensor(cur_fxs);
%        
    cur_dose = doses(i);
    [~,fdose_val] = min(abs(CGobj.mBinsDose - cur_dose));
    vds=zeros(length(cur_grp),1); % volume v at dose d
    vds(:)=0;
    for k=1:length(cur_grp)
        vds(k) = cur_grp(k).fVolAtDose( CGobj.mBinsDose(fdose_val) );
    end

     sa = classKaplanMeierCurve();
  sa.mpValue = inf;
     sa.mHR = 0;
     
     [~,fvol] = min(abs(vds - 31.6));
%     
     flg_volbelow = vds<vds(fvol);
     if sum(flg_volbelow)< 2 || sum(flg_volbelow)>length(vds)-2,
         continue;
     end

     % Cox HR for split
     cox_beta=coxphfit(~flg_volbelow,cur_compdate,'baseline',0,'censoring',cur_flgcensor);
     cox_hr = exp(cox_beta);
      
%     % assign properties of object sa
    survivedate={cur_compdate(flg_volbelow); cur_compdate(~flg_volbelow)}; % survive time of each group
     fcensor={cur_flgcensor(flg_volbelow); cur_flgcensor(~flg_volbelow)}; % censor flag for each group
     sa.mSurvivalTime=survivedate;
     sa.mFlgCensor=fcensor;
%     % compute survival curves and compare them
     sa=sa.fCalculateSurvivalCurve();
     sa=sa.fCombineSurvivalTime();
     sa=sa.fCompareSurvivalByLogrank();
     sa.mHR = cox_hr;
    disp(['Nfx = ',num2str(numFx(i)),' V_{',doses(i),'}']);
    disp(['HR: ',num2str(cox_hr),' p: ',num2str(sa.mpValue)]);
    disp([]);
    
    
      % plot KM curves
    cur_fig=figure(i);  clf reset; hold on; % grid on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'LineWidth',2);
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',12);
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r','LineWidth',2);
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',12);
    ylim([0 0.8]);

    str_pval = ['$p = ',num2str(sa.mpValue,'%10.1e\n'),'$',10,...
            'HR = ',num2str(sa.mHR,3)];
    %text(38,0.25,str_pval,'FontSize',16);
    text(0.25,0.75,str_pval,'FontSize',20,'interpreter','latex');
    lgnd=legend(h_km,...
        ['$V_{',num2str(phys_doses(i)),'} \leq 31.6$cc'],...
        ['$V_{',num2str(phys_doses(i)),'}> 31.6$cc']);
    
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    set(h,'Location','NorthEast');
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Years'],'fontsize',18);
    ylabel(['Probability of CW Pain'],'fontsize',18);
    set(gca,'FontSize',18);
    
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'km_v',num2str(phys_doses(i)),'_fx',num2str(i)],'-pdf');
end
% 
% end

end