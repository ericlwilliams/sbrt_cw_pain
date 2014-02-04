function ChestWallVolumes
tic;
% prepare
fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fn ='MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat';

screen_size=get(0,'ScreenSize');

load(strcat(fp,fn),'CGobj_current');

V30s=zeros(CGobj_current.mNumInGrp,1);
dates=zeros(CGobj_current.mNumInGrp,1);
outcomes=zeros(CGobj_current.mNumInGrp,1);
for k=1:CGobj_current.mNumInGrp
    V30s(k) = CGobj_current.mGrp(k).fVolAtDose( 30 );
    dates(k) = CGobj_current.mGrp(k).mDateStartTx;
    outcomes(k) = CGobj_current.mGrp(k).mFlgCensor;
end


[sorted_dates,idx]=sort(dates);
sorted_v30s = V30s(idx);
sorted_outcomes = outcomes(idx);
sorted_comps = ~sorted_outcomes;

date_range = sorted_dates(end)-sorted_dates(1);
date_width = date_range/39;
binned_dates = sorted_dates+repmat(date_width,[length(sorted_dates) 1]);
binned_comps = histc(sorted_dates(sorted_comps),binned_dates);


% figure(200);
% plot(binned_dates,binned_comps,'o');
% ylabel('Complications/60 days','fontsize',18);
% xlabel('Months','fontsize',18);
% xlim([min(sorted_dates),max(sorted_dates)]);
% datetick('x','mmm/yy','keepticks');
% 


tot_comp = zeros(length(sorted_dates),1);
tot_pts = [1:length(sorted_dates)]';

for i=1:length(sorted_outcomes) %loop over each date
    if i==1
        tot_comp(i)= sorted_comps(i);
    else
        tot_comp(i)=tot_comp(i-1)+sorted_comps(i);
    end
end

% f2=figure(2);
% set(f2,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
% plot(sorted_dates,tot_comp./tot_pts,'LineWidth',2);
% 
% h_ln(1)=line([datenum('07/30/2009') datenum('07/30/2009')],[0 1],...
%     'Color','g',...
%     'LineWidth',2,...
%     'LineStyle','--');
% h_ln(3)=line([datenum('01/01/2011') datenum('01/01/2011')],[0 1],...
%     'Color','b',...
%     'LineWidth',2,...
%     'LineStyle','--');
% h_ln(2)=line([datenum('05/01/2010') datenum('05/01/2010')],[0 1],...
%     'Color','r',...
%     'LineWidth',2,...
%     'LineStyle','--');
% legend(h_ln,'July 30, 2009','May 1, 2010','Jan 1, 2011');
% 
% datetick('x','mmm/yy','keepticks');
% set(gca,'FontSize',12);
% ylabel('Cumulative incidence rate','FontSize',18);
% xlabel('Tx start date','fontsize',18);
% xlim([min(sorted_dates),max(sorted_dates)]);
% ylim([0 1]);


% Get median v30s before/after cut

v30_jan_cut_date='01/01/2011';
v30_jan_cut_date=datenum(v30_jan_cut_date);
v30_bc_jan = sorted_dates<v30_jan_cut_date;

%before
bc_jan_v30s = sorted_v30s(v30_bc_jan);
bc_jan_mdn = median(bc_jan_v30s);

bc_jan_comps = sum(~sorted_outcomes(v30_bc_jan));
bc_jan_outcomes = length(sorted_outcomes(v30_bc_jan));
bc_jan_rate = bc_jan_comps/bc_jan_outcomes;
disp(['Rate of >=2 Gd RP before Jan11: ',num2str(bc_jan_rate,4)]);

% after
ac_jan_v30s = sorted_v30s(~v30_bc_jan);
ac_jan_mdn = median(ac_jan_v30s);

ac_jan_comps = sum(~sorted_outcomes(~v30_bc_jan));
ac_jan_outcomes = length(sorted_outcomes(~v30_bc_jan));
ac_jan_rate = ac_jan_comps/ac_jan_outcomes;
disp(['Rate of >=2 Gd RP after Jan11: ',num2str(ac_jan_rate,4)]);

%may
v30_may_cut_date='05/01/2010';
v30_may_cut_date=datenum(v30_may_cut_date);
v30_bc_may = sorted_dates<v30_may_cut_date;

disp(['']);
%before may
bc_may_v30s = sorted_v30s(v30_bc_may);
bc_may_mdn = median(bc_may_v30s);

bc_may_comps = sum(~sorted_outcomes(v30_bc_may));
bc_may_outcomes = length(sorted_outcomes(v30_bc_may));
bc_may_rate = bc_may_comps/bc_may_outcomes;
disp(['Rate of >=2 Gd RP before May11: ',num2str(bc_may_rate,4)]);

ac_may_v30s = sorted_v30s(~v30_bc_may);
ac_may_mdn = median(ac_may_v30s);

ac_may_comps = sum(~sorted_outcomes(~v30_bc_may));
ac_may_outcomes = length(sorted_outcomes(~v30_bc_may));
ac_may_rate = ac_may_comps/ac_may_outcomes;
disp(['Rate of >=2 Gd RP after May11: ',num2str(ac_may_rate,4)]);



% sept
disp(['Sept 2010']);
v30_sept_cut_date='09/01/2010';
v30_sept_cut_date=datenum(v30_sept_cut_date);
v30_bc_sept = sorted_dates<v30_sept_cut_date;

%before
bc_sept_v30s = sorted_v30s(v30_bc_sept);
bc_sept_mdn = median(bc_sept_v30s);


bc_sept_comps = sum(~sorted_outcomes(v30_bc_sept));
bc_sept_outcomes = length(sorted_outcomes(v30_bc_sept));
bc_sept_rate = bc_sept_comps/bc_sept_outcomes;
disp(['Rate of >=2 Gd RP before Sept10: ',num2str(bc_sept_rate,4)]);

%after

ac_sept_v30s = sorted_v30s(~v30_bc_sept);
ac_sept_mdn = median(ac_sept_v30s);

ac_sept_comps = sum(~sorted_outcomes(~v30_bc_sept));
ac_sept_outcomes = length(sorted_outcomes(~v30_bc_sept));
ac_sept_rate = ac_sept_comps/ac_sept_outcomes;



disp(['Rate of >=2 Gd RP after Sept10: ',num2str(ac_sept_rate,4)]);
disp(['']);
disp(['Median V30 before Sept10: ',num2str(bc_sept_mdn,4)]);
disp(['Median V30 after Sept10: ',num2str(ac_sept_mdn,4)]);


%bc_v30s(bc_v30s==0)=[];

hold on;
f1=figure(1);
set(f1,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
h(1)=plot(sorted_dates,sorted_v30s,'k');
set(gca,'FontSize',12);

xlim([min(sorted_dates),max(sorted_dates)]);
ylabel('Lung volume that received >= 30 Gy [cc]','fontsize',18);
xlabel('Tx start date','fontsize',18);
%set(gca,'XTick',sorted_dates);

h(3)=plot(sorted_dates,repmat(70,length(sorted_dates),1),'g-.','LineWidth',3);

h(4)=plot(sorted_dates(v30_bc_jan),repmat(bc_jan_mdn,length(sorted_dates(v30_bc_jan)),1),'b:','LineWidth',5);
h(5)=plot(sorted_dates(~v30_bc_jan),repmat(ac_jan_mdn,length(sorted_dates(~v30_bc_jan)),1),'b--','LineWidth',3);

h(6)=plot(sorted_dates(v30_bc_may),repmat(bc_may_mdn,length(sorted_dates(v30_bc_may)),1),'r:','LineWidth',5);
h(7)=plot(sorted_dates(~v30_bc_may),repmat(ac_may_mdn,length(sorted_dates(~v30_bc_may)),1),'r--','LineWidth',3);

numBins=78;
binEdges = linspace(min(sorted_dates),max(sorted_dates),numBins+1);



[~,whichBin] = histc(sorted_dates,binEdges);


for i=1:numBins
        flagBinMembers = (whichBin==i);
        binMembersV30 = sorted_v30s(flagBinMembers);
        binMembersDate = sorted_dates(flagBinMembers);
        binMeanV30(i) = mean(binMembersV30);
        binMeanDate(i) = mean(binMembersDate);

end
wideBinEdges = linspace(min(sorted_dates),max(sorted_dates),(numBins/2)+1);
[~,wideWhichBin] = histc(sorted_dates,wideBinEdges);
binNumPts = zeros(numBins/2,1);
for j=1:(numBins/2)
        flagWideBinMembers = (wideWhichBin==j);
        
        wideBinMembersDate = sorted_dates(flagWideBinMembers);
        binOutcomeDate(j)=mean(wideBinMembersDate);
        binOutcomes = sum(sorted_outcomes(flagWideBinMembers))+sum(~sorted_outcomes(flagWideBinMembers));
        binComps = sum(~sorted_outcomes(flagWideBinMembers));
        binCompRate(j) = binComps/binOutcomes;
        if j==1
            binNumPts(j) = binOutcomes;
        else
            binNumPts(j) = binNumPts(j-1)+binOutcomes;
        end
end

binCompRate(isinf(binCompRate))=0;

h(2)=plot(binMeanDate,binMeanV30,'m-','LineWidth',1.5);
yL = get(gca,'YLim');

line([datenum('05/01/2010') datenum('05/01/2010')],yL,...
    'Color','r',...
    'LineStyle','--');
line([datenum('01/01/2011') datenum('01/01/2011')],yL,...
    'Color','b',...
    'LineStyle','--');

datetick('x','mmm/yy','keepticks');

lgnd=legend([h(2) h(3) h(6) h(7) h(4) h(5)], 'Avg./~month','70 cc',...
    strcat('Pre May10 (',num2str(bc_may_mdn,3),' cc)'),...
    strcat('Post May10 (',num2str(ac_may_mdn,3),' cc)'),...
    strcat('Pre Jan11 (',num2str(bc_jan_mdn,3),' cc)'),...    
    strcat('Post Jan11 (',num2str(ac_jan_mdn,3),' cc)'),...
    'Location','NorthEast');
set(lgnd,'FontSize',12);

f100=figure(100);
set(f100,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%plot(binOutcomeDate,binCompRate,'-o','LineWidth',2);
% plot(binOutcomeDate,binNumPts,'-o','LineWidth',2);
% ylabel('Complication Rate/60 days','fontsize',18);
% xlabel('Months','fontsize',18);
% xlim([sorted_dates(1),sorted_dates(end)]);
% set(gca,'XTick',linspace(sorted_dates(1),sorted_dates(end),18));

[AX,H1,H2]=plotyy(binOutcomeDate,binCompRate,binOutcomeDate,binNumPts,'plot');
axes(AX(1));
ylabel('Complication Incidence / 60 days','fontsize',16);
xlim([sorted_dates(1),sorted_dates(end)]);
datetick('x','mmm/yy','keepticks','keeplimits');
axes(AX(2));
ylabel('Cumulative Total Number of Patients','fontsize',16);
xlim([sorted_dates(1),sorted_dates(end)]);
datetick('x','mmm/yy','keepticks','keeplimits');


yL = get(gca,'YLim');
l_jan2011=line([datenum('01/01/2011') datenum('01/01/2011')],yL,...
    'Color','k',...
    'LineWidth',2,...
    'LineStyle','--');
l_jan2010=line([datenum('01/01/2010') datenum('01/01/2010')],yL,...
    'Color','r',...
    'LineWidth',2,...    
    'LineStyle','--');
l_sept2010=line([datenum('09/01/2010') datenum('09/01/2010')],yL,...
    'Color','g',...
    'LineWidth',2,...
    'LineStyle','--');
l_aug2009=line([datenum('07/30/2009') datenum('07/30/2009')],yL,...
    'Color','m',...
    'LineWidth',2,...
    'LineStyle','--');

h_lgnd=legend([l_aug2009 l_jan2010 l_sept2010 l_jan2011],'July 30 2009','Jan 1 2010','Sept 1 2010','Jan 1 2011');
set(h_lgnd,'FontSize',14);

end