function CWP_Plot_a2b

fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
do_print =true;
fig_loc = 'Z:\elw\MATLAB\cw_analy\slides\figures\latest\';


fn=['ChestWall_a2b.mat'];
%fn=['ChestWall_a2b_san.mat'];

%a2b_str = '2.1';
a2b_str = 'Inf';
fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b',a2b_str,'.mat'];
%fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'];
%fn2 = ['RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b0.mat'];

load(strcat(fp,fn));

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

a2b_inf_logl = a2b_logl(1);
a2b_logl = a2b_logl(2:end);

a2b = a2b(2:end);


figure(1); clf reset; hold on;
set(gcf,'Position',ss_four2three);

plot(a2b,a2b_logl,'ko-');hold on;
%[ax,h1,h2]=plotyy(a2b,a2b_logl,a2b,a2b_vx(2:end));hold on;

[max_logl, a2b_max_ind] = max(a2b_logl);
set(gca,'Xtick',0:2:24)

ylabel('Log-likelihood for best V_{NTD} Cox Model','FontSize',20);
xlabel('\alpha/\beta [Gy]','FontSize',20);
%line_sig=plot(xlim,[max_logl-0.5 max_logl-0.5],'r--','LineWidth',2)


set(gca,'FontSize',18)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparing to physical dose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a2b_68cl = interp1(a2b_logl,a2b,a2b_inf_logl+0.5);
if isnan(a2b_68cl)
    a2b_68cl = max(a2b)
end
a2b_95cl = interp1(a2b_logl,a2b,a2b_inf_logl+0.5*chi2inv(0.95,1));

disp(['!!!']);
disp(['Values of a2b that are significantly better V_{D} than PHYSICAL fits at:']);
disp(['68% CI ',num2str(max_logl-0.5),': a2b <= ',num2str(a2b_68cl)]);
disp(['95% CI ',num2str(max_logl-0.5*chi2inv(0.95,1)),': a2b <= ',num2str(a2b_95cl)]);
disp(['!!!']);

%line_68cl=plot([a2b_68cl a2b_68cl],ylim,'g--','LineWidth',2);


%line_mx=plot([a2b(a2b_max_ind) a2b(a2b_max_ind)], ylim,'b--','LineWidth',2);
mx_vxs = vx_logls{a2b_max_ind+1};

[~,vx_ind] = max(mx_vxs);
disp(['LQ Dose for best V_{D_{LQ}} when a2b = ',num2str(a2b(a2b_max_ind))]);
disp(['D_{LQ} = ',num2str(vx_ind)]);

best_dose = vx_ind;

text(9,-317.9,['\alpha/\beta = \infty, best fit LogL = ',...
num2str(a2b_inf_logl,5)],...
     'FontSize',18,...
     'EdgeColor','b',...
     'BackgroundColor','w',...
     'LineWidth',2);


%lgnd=legend([line_68cl line_mx],...
%        ['Upper 68% \alpha/\beta = ',num2str(a2b_68cl,4),' with LogL = ',...
%            num2str(max_logl-0.5,5)],...
%        ['Best fit \alpha/\beta = ',num2str(a2b(a2b_max_ind)),' with LogL = ',...
%            num2str(max_logl,4)],...
%        'Location','Best');
%set(lgnd,'FontSize',18);


 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_vd_vs_a2b'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_vd_vs_a2b.pdf...']);
 end


colors = flipud(varycolor(length(a2b)));

figure(2); clf reset; hold on;
set(gcf,'Position',ss_four2three);
hold on;
for j=1:length(a2b)
    p=plot(1:200,vx_logls{j+1});
    set(p,'Color',colors(j,:));
end
plot(1:200,vx_logls{1},'k','LineWidth',3);
set(gca,'FontSize',18);
ylim([min(vx_logls{2})-2 max(vx_logls{2})+1]);
%grid on;
xlabel('V_{D} [Gy]','FontSize',20)
ylabel('Log-likelihood, Cox model','FontSize',20)
phys_xlim =[xlim];
phys_xlim=phys_xlim(2);
%phys_xlim = xlim(2);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_llhd_vs_vd'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_llhd_vs_vd.pdf...']);
 end
 
 
 %% phys dose only
[best_vd_phys,best_vd_phys_ind] = max(vx_logls{1}); 
 figure(3); clf reset; hold on;
set(gcf,'Position',ss_four2three);
plot(1:200,vx_logls{1},'k','LineWidth',3);hold on;

set(gca,'FontSize',18);
ylim([min(vx_logls{2})-2 max(vx_logls{2})+1]);
h_phys=plot([best_vd_phys_ind best_vd_phys_ind],ylim,'r--','LineWidth',2);
xlim([0 phys_xlim]);
lgnd=legend(h_phys,['Max ln(L) at V_{',num2str(best_vd_phys_ind),'}',10,'is ',num2str(best_vd_phys,4)],'location','best');
set(lgnd,'FontSize',26);
%grid on;
xlabel('V_{D} [Gy]','FontSize',20)
ylabel('Log-likelihood, Cox model','FontSize',20)

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_llhd_vs_vd_phys'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_llhd_vs_vd_phys.pdf...']);
 end
 
 
 figure(100); clf reset; hold on;
set(gcf,'Position',ss_four2three);

plot(a2b,repmat(best_vd_phys,1,length(a2b)),'r--','LineWidth',2);
set(gca,'Xtick',0:2:24)
%ylim([-320.5 -317.5])
ylabel('Log-likelihood for best V_{NTD} Cox Model','FontSize',20);
xlabel('\alpha/\beta [Gy]','FontSize',20);
%line_sig=plot(xlim,[max_logl-0.5 max_logl-0.5],'r--','LineWidth',2)

text(9,-319.5,['Physical Dose',10,'Best fit ln(L) = ',...
num2str(a2b_inf_logl,5)],...
     'FontSize',20,...
     'EdgeColor','b',...
     'BackgroundColor','w',...
     'LineWidth',2);
set(gca,'FontSize',18)

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_phys_llhd_vs_a2b'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_phys_llhd_vs_a2b.pdf...']);
 end
 
 

figure(200); clf reset; hold on;
set(gcf,'Position',ss_four2three);
%plot(a2b,repmat(best_vd_phys,1,length(a2b)),'r--','LineWidth',2);
hold on;

for j=1:length(a2b)
    p=plot(a2b(j),a2b_logl(j),'o');
    set(p,'Color',colors(j,:));
end
[mx_logl,mx_logl_ind]=max(a2b_logl);
low68 = mx_logl-0.5;
h_low68=plot(a2b,repmat(low68,1,length(a2b)),'g--','LineWidth',2);

plot([a2b(mx_logl_ind) a2b(mx_logl_ind)],ylim,'r--','LineWidth',2);

set(gca,'FontSize',18);
set(gca,'Xtick',0:2:24)

text(9,-318.5,['Physical Dose',10,'Best fit ln(L) = ',...
num2str(a2b_inf_logl,5)],...
     'FontSize',20,...
     'EdgeColor','b',...
     'BackgroundColor','w',...
     'LineWidth',2);
 
text(9,-317.9,['NTD Dose',10,'Best fit ln(L) = ',...
num2str(max(a2b_logl),5),10,'at \alpha/\beta = ',num2str(a2b(mx_logl_ind))],...
     'FontSize',20,...
     'EdgeColor','r',...
     'BackgroundColor','w',...
     'LineWidth',2);

lgnd=legend(h_low68,'Low 68% CI','Location','East');
set(lgnd,'FontSize',20);
%ylim([min(vx_logls{2})-2 max(vx_logls{2})+1]);
%grid on;
xlabel('\alpha/\beta [Gy]','FontSize',20)
ylabel('Log-likelihood for best V_{NTD} CPHM','FontSize',20)
%phys_xlim =[xlim];
%phys_xlim=phys_xlim(2);
%phys_xlim = xlim(2);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_llhd_vs_a2b_color'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_llhd_vs_a2b_color.pdf...']);
 end
 
 
 figure(201); clf reset; hold on;
set(gcf,'Position',ss_four2three);
%plot(a2b,repmat(best_vd_phys,1,length(a2b)),'r--','LineWidth',2);
hold on;

[mx_logl,mx_logl_ind]=max(a2b_logl);
best_line_color = [1 0 0];
for j=1:length(a2b)
    p=plot(a2b(j),a2b_logl(j),'o');
    set(p,'Color',colors(j,:));
    if mx_logl_ind==j
        best_line_color = colors(j,:);
    end
end


low68 = mx_logl-0.5*chi2inv(0.68,1);
low95 = mx_logl-0.5*chi2inv(0.95,1);

h_low68=plot(a2b,repmat(low68,1,length(a2b)),'g--','LineWidth',2);
h_low95=plot(a2b,repmat(low95,1,length(a2b)),'b--','LineWidth',2);
h_phys=plot(a2b,repmat(a2b_inf_logl,1,length(a2b)),'k--','LineWidth',3);

% h_best=plot([a2b(mx_logl_ind) a2b(mx_logl_ind)],ylim,'r--','LineWidth',2);
h_best=plot([a2b(mx_logl_ind) a2b(mx_logl_ind)],ylim,'--',...
    'LineWidth',2);
set(h_best,'Color',best_line_color);

set(gca,'FontSize',18);
set(gca,'Xtick',0:2:24)

% text(5,-318.6,['Physical Dose',10,'Best fit ln(L) = ',...
% num2str(a2b_inf_logl,5)],...
%      'FontSize',20,...
%      'EdgeColor','k',...
%      'BackgroundColor','w',...
%      'LineWidth',2);
 
% text(5,-319.1,['NTD Dose',10,'Best fit ln(L) = ',...
% num2str(max(a2b_logl),5),10,'at \alpha/\beta = ',num2str(a2b(mx_logl_ind))],...
%      'FontSize',20,...
%      'EdgeColor','r',...
%      'BackgroundColor','w',...
%      'LineWidth',2);
set(gcf,'Units','normalized')
lgnd=legend([h_best h_low68 h_low95 h_phys],...
    ['Best fit ln(L): ',num2str(max(a2b_logl),5),10,...
    ' at \alpha/\beta = 2.1 Gy'],...
    ['Low 68% CI: ',num2str(max_logl-0.5,5)],...
    ['Low 95% CI: ',num2str(max_logl-0.5*chi2inv(0.95,1),5)],...
    ['Physical dose ln(L): ',num2str(a2b_inf_logl,5)],...
    'Location',[0.3 0.4 0.25 0.25]);

set(lgnd,'FontSize',20);
%ylim([min(vx_logls{2})-2 max(vx_logls{2})+1]);
%grid on;
xlabel('\alpha/\beta [Gy]','FontSize',20)
ylabel('Log-likelihood for best V_{NTD} CPHM','FontSize',20)
%phys_xlim =[xlim];
%phys_xlim=phys_xlim(2);
%phys_xlim = xlim(2);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_llhd_vs_a2b_color_wphys'],'-pdf');
   disp(['Saving ',fig_loc,'cwp_llhd_vs_a2b_color_wphys.pdf...']);
 end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
%% plot KM curve for best a2b (= 2.1)

return;

load(strcat(fp,fn2),'CGobj_current');
CGobj = {CGobj_current};


split = -1;
dose=best_dose;    
%dose=42;    
[cur_fig,~,~] = fPlotKaplanMeierCurve_VDx(CGobj,dose,split);
 if do_print,
    set(cur_fig,'Color','w');
    export_fig(gcf,[fig_loc,'cwp_km_med_split_a2b',a2b_str],'-pdf');
   disp(['Saving ',fig_loc,'cwp_km_med_split_a2b',a2b_str,'.pdf...']);
 end

split = -1;
dose=best_dose;    
fPlotKaplanMeierCurve(CGobj,dose,split,3);

split = -1;
dose=best_dose;    
fPlotKaplanMeierCurve(CGobj,dose,split,4);

split = -1;
dose=best_dose;
fPlotKaplanMeierCurve(CGobj,dose,split,5);

%best split at 108cc (within 2-3rd quartile)
%(determined in ChestWallLogRankDisplay.m (have to change config)
split = 108;
dose=best_dose;
fPlotKaplanMeierCurve(CGobj,dose,split);


%best split mid quartiles only is: 
% split = 108;
% dose=99;    
% fPlotKaplanMeierCurve(CGobj,dose,split);

end
