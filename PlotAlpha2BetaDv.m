function PlotAlpha2BetaDv

%fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fig_loc = 'Z:\elw\MATLAB\cw_analy\slides\figures\latest\';

fn=['ChestWall_a2b_dv.mat'];
%fn=['ChestWall_a2b_san.mat'];
do_print=1;

fn2 = ['MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b2.1.mat'];
%fn2 = ['RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'];
%fn2 = ['RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2b0.mat'];

load(strcat(fp,fn));

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

a2b_inf_logl = a2b_logl(1);
a2b_logl = a2b_logl(2:end);

a2b = a2b(2:end);


cur_fig=figure(1); clf reset; hold on;
set(gcf,'Position',ss_four2three);

plot(a2b,a2b_logl,'ko-');hold on;
%[ax,h1,h2]=plotyy(a2b,a2b_logl,a2b,a2b_dx(2:end));hold on;

[max_logl, a2b_max_ind] = max(a2b_logl);

ylabel('Log-likelihood for best D_{V} Cox Model','FontSize',14);
xlabel('\alpha/\beta [Gy]','FontSize',14);
line_sig=plot(xlim,[max_logl-0.5 max_logl-0.5],'r--','LineWidth',2)
%line_sig=plot(xlim,[a2b_inf_logl+1.96 a2b_inf_logl+1.96],'r--','LineWidth',2)
grid on;


set(gca,'FontSize',12)

%line_inf = plot(xlim,[a2b_inf_logl a2b_inf_logl],'m--','LineWidth',2);

% below a2b_ul, fit is significantly improved than w/out a2b corr

a2b_68cl = interp1(a2b_logl,a2b,max_logl-0.5);

a2b_ul = interp1(a2b_logl,a2b,a2b_inf_logl+1.96);
disp(['a2b <= ',num2str(a2b_ul),' significantly better V_{D} fit than no a2b']);
%line_ul=plot([a2b_ul a2b_ul],ylim,'g--','LineWidth',2);
line_68cl=plot([a2b_68cl a2b_68cl],ylim,'g--','LineWidth',2);


line_mx=plot([a2b(a2b_max_ind) a2b(a2b_max_ind)], ylim,'b--','LineWidth',2);
mx_dxs = dx_logls{a2b_max_ind+1};

[~,dx_ind] = max(mx_dxs);
disp(['LQ Dose for best V_{D_{LQ}} when a2b = ',num2str(a2b(a2b_max_ind))]);
disp(['D_{LQ} = ',num2str(dx_ind)]);

best_dose = dx_ind;

%text(9,-317.9,['\alpha/\beta = \infty, best fit LogL = ',...
text(7,-322.5,['\alpha/\beta = \infty, best fit LogL = ',...
    num2str(a2b_inf_logl,5)],...
     'FontSize',16,...
     'EdgeColor','b',...
     'BackgroundColor','w',...
     'LineWidth',2);

%lgnd=legend([line_sig line_68cl line_mx],...
        %['Lower 68% CL for best fit LogL = ',num2str(max_logl-0.5,5)],...
lgnd=legend([line_68cl line_mx],...
        ['Upper 68% \alpha/\beta = ',num2str(a2b_68cl,4),' with LogL = ',...
            num2str(max_logl-0.5,5)],...
        ['Best fit \alpha/\beta = ',num2str(a2b(a2b_max_ind)),' with LogL = ',...
            num2str(max_logl,4)],...
        'Location','Best');
set(lgnd,'FontSize',14);



 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'dv_logl_vs_a2b'],'-pdf');
   disp(['Saving ',fig_loc,'dv_logl_vs_a2b.pdf...']);
 end



cur_fig=figure(2); clf reset; hold on;
set(gcf,'Position',ss_four2three);

colors = flipud(varycolor(length(a2b)));

hold on;
for j=1:length(a2b)
    p=plot(1:200,dx_logls{j+1});
    set(p,'Color',colors(j,:));
end
plot(1:200,dx_logls{1},'k','LineWidth',3);

ylim([min(dx_logls{2})-2 max(dx_logls{2})+1]);
grid on;
xlabel('D_{V} [cc]','FontSize',14)
ylabel('Log-likelihood, Cox model','FontSize',14)

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'logl_vs_dv'],'-pdf');
   disp(['Saving ',fig_loc,'logl_vs_dv.pdf...']);
 end
%% plot KM curve for best a2b (= 2.1)

return;
% 
% load(strcat(fp,fn2),'CGobj_current');
% CGobj = {CGobj_current};
% 
% 
% split = -1;
% dose=best_dose;    
% %dose=42;    
% fPlotKaplanMeierCurve(CGobj,dose,split);
% 
% split = -1;
% dose=best_dose;    
% fPlotKaplanMeierCurve(CGobj,dose,split,3);
% 
% split = -1;
% dose=best_dose;    
% fPlotKaplanMeierCurve(CGobj,dose,split,4);
% 
% split = -1;
% dose=best_dose;
% fPlotKaplanMeierCurve(CGobj,dose,split,5);
% 
% %best split at 108cc (within 2-3rd quartile)
% %(determined in ChestWallLogRankDisplay.m (have to change config)
% split = 108;
% dose=best_dose;
% fPlotKaplanMeierCurve(CGobj,dose,split);
% 
% 
% %best split mid quartiles only is: 
% % split = 108;
% % dose=99;    
% % fPlotKaplanMeierCurve(CGobj,dose,split);

end
