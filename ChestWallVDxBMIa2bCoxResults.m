function ChestWallVDxBMIa2bCoxResults
tic;

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fn = 'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat';
screen_size=get(0,'ScreenSize');


% load data
load(strcat(fp,fn),'CGobj_current');
CGobj = CGobj_current;

% select patients with data
f = CGobj.fPatientsWithComplicationData();
CGobj = CGobj.fRemovePatient(~f);

f = cellfun(@(x) strcmpi('VDxBMI',x),CGobj.mCoxParameter(:,1)); % search the label
VDxBMICox = CGobj.mCoxParameter{f,2}; % extract Cox model result
flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), VDxBMICox); % some fields are empty or infinite, indicating no data for those values
VDxBMICox = VDxBMICox(flgCox);
tmp_beta = {VDxBMICox.beta};
is_inf = cellfun(@(x) length(x),tmp_beta);
is_inf = is_inf==1;
tmp_beta = tmp_beta(~is_inf);


VDxBMICox = VDxBMICox(~is_inf);
betas = cell2mat(tmp_beta)';


v99_idx=100;
VD99BMICox = VDxBMICox(v99_idx);
cur_vd = VD99BMICox.data_exposure(:,1);
cur_vd_beta = betas(v99_idx,1);

cur_bmi = VD99BMICox.data_exposure(:,2);
cur_bmi_beta = betas(v99_idx,2);

cur_arg = [cur_vd_beta*cur_vd+cur_bmi_beta*cur_bmi];

splits = cur_arg;

hzrd = VD99BMICox.h;
compdate = VD99BMICox.data_hazard;
flgcensor = [CGobj.mGrp.mFlgCensor]';
bmi = [CGobj.mGrp.mBMI]';
bmi_idx = bmi>0;
flgcensor = flgcensor(bmi_idx);


 f = find(diff(hzrd(:,1))==0); % find duplicate time values of h(t)
 while ~isempty(f)
     hzrd(f,1) = hzrd(f,1)-eps*10; % adjust it a bit to avoid ambiguius
     f = find(diff(hzrd(:,1))==0); % find duplicate time values of h(t)
 end
 
fig=figure(1); clf reset;hold on;
set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
 

grid on;
asympt_rates = zeros(length(splits),1);
raw_rates = zeros(length(splits),1);
cox_z = zeros(length(splits),1);

[splits,split_idx] = sort(splits); 
compdate = compdate(split_idx); 
flgcensor = flgcensor(split_idx); 
%splits = unique(splits);

colors = varycolor(length(splits));

flg_split2=-1;
numstart=4;
numend=length(splits)-numstart;
skip_split = zeros(length(splits),1);
cur_color=1;
for i=1:length(splits)

    cur_split = splits(i);  
    
    %flg_split=cur_split<=splits;
    flg_split=splits<=cur_split;
    f=length(find(flg_split));
    if f<numstart || f>numend % one group has too few patients, skip it
        skip_split(i)=1;
        continue;
    end
    if isequal(flg_split,flg_split2) % if it is the same grouping as previous, skip the computation and save the result directly
        skip_split(i)=1;
        continue;
    end
    
    flg_split2=flg_split;
    
    [cur_b,~,cur_H,cur_stats] =...
        coxphfit(flg_split,compdate,'baseline',0,'censoring',flgcensor);

    if cur_stats.p>0.05
        skip_split(i)=1;
        continue;
    end
    
    expbetax_above = 1;
    expbetax_below = exp(cur_b);
    cur_cic_above = 1-exp(-cur_H(:,2)*expbetax_above);
    cur_cic_below = 1-exp(-cur_H(:,2)*expbetax_below);

    if cur_cic_below(end) > cur_cic_above(end)
        skip_split(i)=1;
        continue
    end
    
   p=plot(cur_H(:,1), cur_cic_below,'b-',...
        'LineWidth',1); % high resolution Cox curve
   set(p,'Color',colors(cur_color,:));

   cur_color=cur_color+1;
   
   asympt_rates(i)=cur_cic_below(end);    
    %raw_rates(i)=mean(cur_cic_below);
    raw_rates(i)=...
        sum(flg_split(~flgcensor))/f; %number complications below split/total below split
        
    cox_z(i)=abs(cur_stats.z);
    %cox_z(i)=abs(cur_stats.p);
    

end
 splits(skip_split==1)=[];
 asympt_rates(skip_split==1)=[];
 raw_rates(skip_split==1)=[];
 cox_z(skip_split==1)=[];
 
 
 hold off;
 title('Lower Cumulative CWP Incidence','FontSize',15);
 set(gca,'xminortick','on','yminortick','on');
 set(gca,'box','on');
xlabel('Month','FontSize',15);
ylabel('Probability of Complication','FontSize',15);

  
  %% Lower rates vs split value
fig=figure(2); clf reset;hold on;
set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
grid on;
%h_asympt=plot(splits,asympt_rates,'b-.');
[ax,h_rates,h_z] =plotyy([splits splits],[asympt_rates raw_rates],splits,cox_z);
set(get(ax(1),'Ylabel'),'String','Lower CWP Rate','FontSize',15);
%set(ax(1),'ylim',[0 0.3]);
set(get(ax(2),'Ylabel'),'String','Cox Z-value','FontSize',15);
%set(ax(2),'ylim',[0 7]);

set(h_rates,'LineWidth',2)
set(h_z,'LineWidth',2)
h_lgnd=legend([h_rates;h_z],...
            'Asymptotic Incidence',...
            'Raw Incidence',...
            'Cox Z');
xlabel('Split value','FontSize',15);
%ylabel('Lower CWP Rate','FontSize',15);
set(h_lgnd,'fontsize',14);
set(h_lgnd,'location','best');

 %% Vd vs BMI
b_bmi = cur_bmi_beta;
b_vd = cur_vd_beta;

fig=figure(5); clf reset;hold on;
set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

for j=1:length(splits)
%     if skip_split(j)
%         continue;
%     end
    bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + splits(j)/b_vd;
    h_mid_pct=ezplot(bmi_vs_vd,[0,40]);
    set(h_mid_pct,'Color',colors(j,:));
end
title('Threshold for incidence in lower risk group','FontSize',15);
xlabel(['BMI'],'FontSize',15);
ylabel(['V_{99,a/b=2.1} [cc]'],'FontSize',15); 
ylim([0 225]);
grid on;
hold off;



%% Vd vs BMI with splits grouped by rate, and weigthed by z
% lower  R <= 0.125
% mid    0.125 > R <= 0.175
% high   R > 0.175
low_raw_idx = raw_rates<=quantile(raw_rates,1/3);
mid_raw_idx = [raw_rates>quantile(raw_rates,1/3)].*...
            [raw_rates<=quantile(raw_rates,2/3)];
high_raw_idx = raw_rates>quantile(raw_rates,2/3);

low_splits = splits(low_raw_idx==1);
se_low_splits = std(low_splits)/sqrt(length(low_splits));
low_splits=low_splits.*cox_z(low_raw_idx==1);
mean_low_splits=sum(low_splits)/sum(cox_z(low_raw_idx==1));
low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + mean_low_splits/b_vd;
low_low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_low_splits-se_low_splits)/b_vd;
low_high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_low_splits+se_low_splits)/b_vd;

%mid_idx = [asympt_rates>0.125].*[asympt_rates<=0.175];


mid_splits = splits(mid_raw_idx==1);
se_mid_splits = std(mid_splits)/(sqrt(length(mid_splits)));
mid_splits=mid_splits.*cox_z(mid_raw_idx==1);
mean_mid_splits=sum(mid_splits)/sum(cox_z(mid_raw_idx==1));
mid_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + mean_mid_splits/b_vd;
mid_low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_mid_splits-se_mid_splits)/b_vd;
mid_high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_mid_splits+se_mid_splits)/b_vd;


%high_idx = asympt_rates>0.175;

    
high_splits = splits(high_raw_idx==1);
se_high_splits = std(high_splits)/sqrt(length(high_splits));
high_splits=high_splits.*cox_z(high_raw_idx==1);
mean_high_splits=sum(high_splits)/sum(cox_z(high_raw_idx==1));
high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + mean_high_splits/b_vd;
high_low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_high_splits-se_high_splits)/b_vd;
high_high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
    (mean_high_splits+se_high_splits)/b_vd;


fig=figure(6); clf reset;hold on;
set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_low_pct=ezplot(low_bmi_vs_vd,[0,40]);
set(h_low_pct,'Color','g');
set(h_low_pct,'LineWidth',2);
h_low_low_pct=ezplot(low_low_bmi_vs_vd,[0,40]);
set(h_low_low_pct,'Color','g');
set(h_low_low_pct,'LineStyle',':');
h_low_high_pct=ezplot(low_high_bmi_vs_vd,[0,40]);
set(h_low_high_pct,'Color','g');
set(h_low_high_pct,'LineStyle',':');

h_mid_pct=ezplot(mid_bmi_vs_vd,[0,40]);
set(h_mid_pct,'Color','b');
set(h_mid_pct,'LineWidth',2);
h_mid_low_pct=ezplot(mid_low_bmi_vs_vd,[0,40]);
set(h_mid_low_pct,'Color','b');
set(h_mid_low_pct,'LineStyle',':');
h_mid_high_pct=ezplot(mid_high_bmi_vs_vd,[0,40]);
set(h_mid_high_pct,'Color','b');
set(h_mid_high_pct,'LineStyle',':');

h_high_pct=ezplot(high_bmi_vs_vd,[0,40]);
set(h_high_pct,'Color','r');
set(h_high_pct,'LineWidth',2);
h_high_low_pct=ezplot(high_low_bmi_vs_vd,[0,40]);
set(h_high_low_pct,'Color','r');
set(h_high_low_pct,'LineStyle',':');
h_high_high_pct=ezplot(high_high_bmi_vs_vd,[0,40]);
set(h_high_high_pct,'Color','r');
set(h_high_high_pct,'LineStyle',':');


ylim([0 160]);
%ylim([0 210]);
title(['Threshold for raw incidence in lower risk group',10,...
    'Weighted by Cox z-values, errors \pm 1 S.E.'],'FontSize',15);
xlabel(['BMI'],'FontSize',15);
ylabel(['V_{99,a/b=2.1} [cc]'],'FontSize',15); 
grid on;

lgnd=legend([h_low_pct h_mid_pct h_high_pct],...
        ['R_{raw} \leq ',num2str(quantile(raw_rates,1/3)*100,3),'%'],...
        [num2str(quantile(raw_rates,1/3)*100,3),'% < R_{raw} \leq ',...
        num2str(quantile(raw_rates,2/3)*100,3),'%'],...
        ['R_{raw} > ',num2str(quantile(raw_rates,2/3)*100,3),'%']);

set(lgnd,'fontsize',14);
set(lgnd,'location','best');
 %% vd vs bmi

fig=figure(3); clf reset;hold on;
%set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
[n_low,x_low]=hist(raw_rates(low_raw_idx==1),0:0.01:0.3);
bar(x_low,n_low,'g');

[n_mid,x_mid]=hist(raw_rates(mid_raw_idx==1),0:0.01:0.3);
bar(x_mid,n_mid,'b');

[n_high,x_high]=hist(raw_rates(high_raw_idx==1),0:0.01:0.3);
bar(x_high,n_high,'r');

xlim([0 0.32]);
xlabel('Lower raw CWP incidence','FontSize',15);
ylabel('Frequency','FontSize',15);


% % asympt
% fig=figure(7); clf reset;hold on;
% %set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
% [n_low,x_low]=hist(asympt_rates(low_idx==1),0:0.01:0.3);
% bar(x_low,n_low,'g');
% 
% [n_mid,x_mid]=hist(asympt_rates(mid_idx==1),0:0.01:0.3);
% bar(x_mid,n_mid,'b');
% 
% [n_high,x_high]=hist(asympt_rates(high_idx==1),0:0.01:0.3);
% bar(x_high,n_high,'r');
% 
% xlim([0 0.32]);
% xlabel('Lower asymptotic CWP incidence','FontSize',15);
% ylabel('Frequency','FontSize',15);



%% See raw < 10 % and asympt < 10%
% find best split for 9%<R < 11$
raw_10pct_idx = [raw_rates>0.09].*[raw_rates<0.11];
%[raw_mx_z_10pct,~] = max(cox_z(raw_10pct_idx==1));
%raw_mx_z_idx =  find(cox_z==raw_mx_z_10pct);
%raw_10pct_split = splits(raw_mx_z_idx);
raw_10pct_split = max(splits(raw_10pct_idx==1));
raw_10pct_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + raw_10pct_split/b_vd;

asympt_10pct_idx = [asympt_rates>0.09].*[asympt_rates<0.11];
%[asympt_mx_z_10pct,~] = max(cox_z(asympt_10pct_idx==1));
%asympt_mx_z_idx =  find(cox_z==asympt_mx_z_10pct);
%asympt_10pct_split = splits(asympt_mx_z_idx);

asympt_10pct_split = max(splits(asympt_10pct_idx==1));
asympt_10pct_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + asympt_10pct_split/b_vd;


fig=figure(8); clf reset;hold on;
set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_raw_10pct=ezplot(raw_10pct_bmi_vs_vd,[0,40]);
set(h_raw_10pct,'Color','r');
set(h_raw_10pct,'LineWidth',2);
h_asympt_10pct=ezplot(asympt_10pct_bmi_vs_vd,[0,40]);
set(h_asympt_10pct,'Color','b');
set(h_asympt_10pct,'LineWidth',2);

ylim([0 125]);
title('Max V_{99} threshold for 10% incidence in lower risk group',...
    'FontSize',15);
xlabel(['BMI'],'FontSize',15);
ylabel(['V_{99,a/b=2.1} [cc]'],'FontSize',15); 
grid on;

lgnd=legend([h_raw_10pct h_asympt_10pct],...
        ['R_{raw} \leq 10%'],...
        ['R_{asympt} \leq 10%']);
set(lgnd,'fontsize',14);
set(lgnd,'location','best');


% 
% raw_10pct_idx = raw_rates<=0.1;
% asympt_10pct_idx = asympt_rates<=0.1;
% 
% 
% raw_10pct_splits = splits(raw_10pct_idx==1);
% se_raw_10pct_splits = std(raw_10pct_splits)/sqrt(length(raw_10pct_splits));
% raw_10pct_splits=raw_10pct_splits.*cox_z(raw_10pct_idx==1);
% mean_raw_10pct_splits=sum(raw_10pct_splits)/sum(cox_z(raw_10pct_idx==1));
% 
% raw_10pct_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + mean_raw_10pct_splits/b_vd;
% raw_10pct_low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
%     (mean_raw_10pct_splits-se_raw_10pct_splits)/b_vd;
% raw_10pct_high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
%     (mean_raw_10pct_splits+se_raw_10pct_splits)/b_vd;
% 
% asympt_10pct_splits = splits(asympt_10pct_idx==1);
% se_asympt_10pct_splits = std(asympt_10pct_splits)/sqrt(length(asympt_10pct_splits));
% asympt_10pct_splits=asympt_10pct_splits.*cox_z(asympt_10pct_idx==1);
% mean_asympt_10pct_splits=sum(asympt_10pct_splits)/sum(cox_z(asympt_10pct_idx==1));
% 
% asympt_10pct_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x + mean_asympt_10pct_splits/b_vd;
% asympt_10pct_low_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
%     (mean_asympt_10pct_splits-se_asympt_10pct_splits)/b_vd;
% asympt_10pct_high_bmi_vs_vd = @(x) (-b_bmi/b_vd)*x +...
%     (mean_asympt_10pct_splits+se_asympt_10pct_splits)/b_vd;
% 
% 
% fig=figure(8); clf reset;hold on;
% set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
% 
% 
% 
% h_raw_10pct=ezplot(raw_10pct_bmi_vs_vd,[0,40]);
% set(h_raw_10pct,'Color','r');
% set(h_raw_10pct,'LineWidth',2);
% h_raw_10pct_low=ezplot(raw_10pct_low_bmi_vs_vd,[0,40]);
% set(h_raw_10pct_low,'Color','r');
% set(h_raw_10pct_low,'LineStyle',':');
% h_raw_10pct_high=ezplot(raw_10pct_high_bmi_vs_vd,[0,40]);
% set(h_raw_10pct_high,'Color','r');
% set(h_raw_10pct_high,'LineStyle',':');
% 
% h_asympt_10pct=ezplot(asympt_10pct_bmi_vs_vd,[0,40]);
% set(h_asympt_10pct,'Color','b');
% set(h_asympt_10pct,'LineWidth',2);
% h_asympt_10pct_low=ezplot(asympt_10pct_low_bmi_vs_vd,[0,40]);
% set(h_asympt_10pct_low,'Color','b');
% set(h_asympt_10pct_low,'LineStyle',':');
% h_asympt_10pct_high=ezplot(asympt_10pct_high_bmi_vs_vd,[0,40]);
% set(h_asympt_10pct_high,'Color','b');
% set(h_asympt_10pct_high,'LineStyle',':');
% 
% 
% ylim([0 95]);
% title(['Threshold for 10%raw/asymptotic incidence in lower risk group',10,...
%     'Weighted by Cox z-values, errors \pm 1 S.E.'],'FontSize',15);
% xlabel(['BMI'],'FontSize',15);
% ylabel(['V_{99,a/b=2.1} [cc]'],'FontSize',15); 
% grid on;
% 
% lgnd=legend([h_raw_10pct h_asympt_10pct],...
%         ['R_{raw} \leq 10%'],...
%         ['R_{asympt} \leq 10%']);
% set(lgnd,'fontsize',14);
% set(lgnd,'location','best');



 end