function ChestWallPainAlpha2BetaOLD
tic;

%fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
fp = 'C:\Documents and Settings\williae1\cw_meta_data\';
fn = 'MUTTER_MASTER_ChestWall_CoxPHM_Results_a2b_b200.mat';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
ss_full = screen_size;
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];

%% load data
% produced from WriteCoxModelResults.m
load(strcat(fp,fn));

%% Print DVx logl

a2bInf_logl = [a2bInf_DVxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
[a2bInf_mx,doseloc]=max(a2bInf_logl); % the best fitting of Cox model
a2bInf_lowCI68 = a2bInf_mx - 0.5; % 68% confidence
a2bInf_lowCI95 = a2bInf_mx - 1.96; % 95% confidence

a2b3_logl = [a2b3_DVxCox.logl]'; %logl(~flgCox) = -3; % log likelihood of Cox model, anti-correlation points not counted
[a2b3_mx,doseloc]=max(a2b3_logl); % the best fitting of Cox model
a2b3_lowCI68 = a2b3_mx - 0.5; % 68% confidence
a2b3_lowCI95 = a2b3_mx - 1.96; % 95% confidence


figure(1); clf reset; hold on;
set(gcf,'Position',ss_four2three);
h_a2b3=plot(x_a2b3_dvx, [a2b3_DVxCox.logl],'r.-');
h_a2b15=plot(x_a2b15_dvx, [a2b15_DVxCox.logl],'.-');
h_a2b10=plot(x_a2b10_dvx, [a2b10_DVxCox.logl],'c.-');
h_a2b200=plot(x_a2b200_dvx, [a2b200_DVxCox.logl],'m.-');
h_a2b75=plot(x_a2b75_dvx, [a2b75_DVxCox.logl],'g.-');
h_a2bInf=plot(x_a2bInf_dvx, [a2bInf_DVxCox.logl],'k.-');
xlim([0 max(x_a2b3_dvx)]);
ylim([-338 -321.5]);

plot(x_a2b3_dvx,repmat(a2b3_lowCI68,size(x_a2b3_dvx)),'r--','LineWidth',2);
plot(x_a2b3_dvx,repmat(a2b3_lowCI95,size(x_a2b3_dvx)),'c--','LineWidth',2);
h_a2b3_mx_logl=plot([x_a2b3_dvx(doseloc) x_a2b3_dvx(doseloc)],ylim,'g--','LineWidth',2);



set(gca,'fontsize',14);
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('log likelihood','fontsize',18);
%logl_str = {['Max LogL = ', num2str(a2b3_mx,5),' at ',num2str(x_a2b3_dvx(doseloc),4),' cc']};
%legend(h_a2b3_mx_logl,logl_str);
h_lgnd=legend([h_a2b3 h_a2b10 h_a2b15 h_a2b75 h_a2b200 h_a2bInf],...
    '$\frac{\alpha}{\beta}=3~$Gy',...
    '$\frac{\alpha}{\beta}=10~$Gy',...
    '$\frac{\alpha}{\beta}=15~$Gy',...
    '$\frac{\alpha}{\beta}=75~$Gy',...
    '$\frac{\alpha}{\beta}=200~$Gy',...
    '$\frac{\alpha}{\beta}=~$Inf',...
    'Location','Best');
set(h_lgnd,'fontsize',18);
set(h_lgnd,'interpreter','latex');

%% Print DVx pvals
a2b3_DVxCox_pvals = ones(size(a2b3_DVxCox));
a2b3_DVxCox_pvals = [a2b3_DVxCox.p]';
[min_a2b3_DVxCox_pvals,ploc] = min(a2b3_DVxCox_pvals);

figure(2); clf reset;hold off;
set(gcf,'Position',ss_four2three);
h_a2b3=semilogy(x_a2b3_dvx,[a2b3_DVxCox.p],'r.-');
hold on;
h_a2bInf=semilogy(x_a2bInf_dvx,[a2bInf_DVxCox.p],'k.-');
h_a2b15=semilogy(x_a2b15_dvx,[a2b15_DVxCox.p],'.-');
h_a2b10=semilogy(x_a2b10_dvx,[a2b10_DVxCox.p],'c.-');
h_a2b200=semilogy(x_a2b200_dvx,[a2b200_DVxCox.p],'m.-');
h_a2b75=semilogy(x_a2b75_dvx,[a2b75_DVxCox.p],'g.-');
semilogy(x_a2b75_dvx,repmat(0.05,size(x_a2b75_dvx)),'c--');
xlim([0 1000]);
h_min_pval = semilogy([x_a2b3_dvx(ploc) x_a2b3_dvx(ploc)],ylim,'g--','LineWidth',2);

hold off;
xlim([0 max(x_a2b3_dvx)]);
set(gca,'fontsize',14);
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
pval_str = {['Min p-val = ', num2str(min_a2b3_DVxCox_pvals,3),' at ',...
    num2str(x_a2b3_dvx(ploc),4),' cc']};
h_lgnd=legend([h_a2b3 h_a2b10 h_a2b15 h_a2b75 h_a2b200 h_a2bInf],...
    '$\frac{\alpha}{\beta}=3~$Gy',...
    '$\frac{\alpha}{\beta}=10~$Gy',...
    '$\frac{\alpha}{\beta}=15~$Gy',...
    '$\frac{\alpha}{\beta}=75~$Gy',...
    '$\frac{\alpha}{\beta}=200~$Gy',...
    '$\frac{\alpha}{\beta}=~$Inf',...
    'Location','Best');
set(h_lgnd,'fontsize',18);
set(h_lgnd,'interpreter','latex');

xlabel('(D_{V}) Volume [cc]','fontsize',18); ylabel('p-value','fontsize',18);

%% Print VDx logl

a2bInf_logl = [a2bInf_VDxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
[a2bInf_mx,doseloc]=max(a2bInf_logl); % the best fitting of Cox model
a2bInf_lowCI68 = a2bInf_mx - 0.5; % 68% confidence
a2bInf_lowCI95 = a2bInf_mx - 1.96; % 95% confidence

a2b3_logl = [a2b3_VDxCox.logl]'; %logl(~flgCox) = -3; % log likelihood of Cox model, anti-correlation points not counted
[a2b3_mx,doseloc]=max(a2b3_logl); % the best fitting of Cox model
a2b3_lowCI68 = a2b3_mx - 0.5; % 68% confidence
a2b3_lowCI95 = a2b3_mx - 1.96; % 95% confidence


figure(3); clf reset; hold on;
set(gcf,'Position',ss_four2three);
h_a2b3=plot(x_a2b3_vdx, [a2b3_VDxCox.logl],'r.-');
h_a2b15=plot(x_a2b15_vdx, [a2b15_VDxCox.logl],'.-');
h_a2b10=plot(x_a2b10_vdx, [a2b10_VDxCox.logl],'c.-');
h_a2b200=plot(x_a2b200_vdx, [a2b200_VDxCox.logl],'m.-');
h_a2b75=plot(x_a2b75_vdx, [a2b75_VDxCox.logl],'g.-');
h_a2bInf=plot(x_a2bInf_vdx, [a2bInf_VDxCox.logl],'k.-');
%xlim([0 max(x_a2b3_vdx)]);
xlim([0 150]);
ylim([-337 -317]);

plot(x_a2b3_vdx,repmat(a2b3_lowCI68,size(x_a2b3_vdx)),'r--','LineWidth',2);
plot(x_a2b3_vdx,repmat(a2b3_lowCI95,size(x_a2b3_vdx)),'c--','LineWidth',2);
h_a2b3_mx_logl=plot([x_a2b3_vdx(doseloc) x_a2b3_vdx(doseloc)],ylim,'g--','LineWidth',2);
hold off;


set(gca,'fontsize',14);
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
xlabel('(V_{D}) Dose [Gy]','fontsize',18); ylabel('log likelihood','fontsize',18);
%logl_str = {['Max LogL = ', num2str(a2b3_mx,5),' at ',num2str(x_a2b3_vdx(doseloc),4),' cc']};
h_lgnd=legend([h_a2b3 h_a2b10 h_a2b15 h_a2b75 h_a2b200 h_a2bInf],...
    '$\frac{\alpha}{\beta}=3~$Gy',...
    '$\frac{\alpha}{\beta}=10~$Gy',...
    '$\frac{\alpha}{\beta}=15~$Gy',...
    '$\frac{\alpha}{\beta}=75~$Gy',...
    '$\frac{\alpha}{\beta}=200~$Gy',...
    '$\frac{\alpha}{\beta}=~$Inf',...
    'Location','Best');
set(h_lgnd,'fontsize',18);
set(h_lgnd,'interpreter','latex');

%here
a2b15_logl = [a2b15_VDxCox.logl]'; %logl(~flgCox) = -3; % log likelihood of Cox model, anti-correlation points not counted
[a2b15_mx,doseloc]=max(a2b15_logl); % the best fitting of Cox model
a2b15_lowCI68 = a2b15_mx - 0.5; % 68% confidence
a2b15_lowCI95 = a2b15_mx - 1.96; % 95% confidence


figure(30); clf reset; hold on;
set(gcf,'Position',ss_four2three);
h_a2b3=plot(x_a2b3_vdx, [a2b3_VDxCox.logl],'r.-','LineWidth',2);
h_a2b15=plot(x_a2b15_vdx, [a2b15_VDxCox.logl],'.-','LineWidth',2);
h_a2b10=plot(x_a2b10_vdx, [a2b10_VDxCox.logl],'c.-','LineWidth',2);
h_a2b200=plot(x_a2b200_vdx, [a2b200_VDxCox.logl],'m.-','LineWidth',2);
h_a2b75=plot(x_a2b75_vdx, [a2b75_VDxCox.logl],'g.-','LineWidth',2);
h_a2bInf=plot(x_a2bInf_vdx, [a2bInf_VDxCox.logl],'k.-','LineWidth',2);
%xlim([0 max(x_a2b3_vdx)]);
xlim([0 115]);
ylim([-322 -317.5]);

plot(x_a2b15_vdx,repmat(a2b15_lowCI68,size(x_a2b15_vdx)),'r--','LineWidth',2);
plot(x_a2b15_vdx,repmat(a2b15_lowCI95,size(x_a2b15_vdx)),'c--','LineWidth',2);
h_a2b15_mx_logl=plot([x_a2b15_vdx(doseloc) x_a2b15_vdx(doseloc)],ylim,'g--','LineWidth',2);
hold off;


set(gca,'fontsize',14);
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
xlabel('(V_{D}) Dose [Gy]','fontsize',18); ylabel('log likelihood','fontsize',18);
%logl_str = {['Max LogL = ', num2str(a2b3_mx,5),' at ',num2str(x_a2b3_vdx(doseloc),4),' cc']};
h_lgnd=legend([h_a2b3 h_a2b10 h_a2b15 h_a2b75 h_a2b200 h_a2bInf],...
    '$\frac{\alpha}{\beta}=3~$Gy',...
    '$\frac{\alpha}{\beta}=10~$Gy',...
    '$\frac{\alpha}{\beta}=15~$Gy',...
    '$\frac{\alpha}{\beta}=75~$Gy',...
    '$\frac{\alpha}{\beta}=200~$Gy',...
    '$\frac{\alpha}{\beta}=~$Inf',...
    'Location','Best');
set(h_lgnd,'fontsize',18);
set(h_lgnd,'interpreter','latex');

%% Print VDx pvals
a2b3_VDxCox_pvals = ones(size(a2b3_VDxCox));
a2b3_VDxCox_pvals = [a2b3_VDxCox.p]';
[min_a2b3_VDxCox_pvals,ploc] = min(a2b3_VDxCox_pvals);

figure(4); clf reset;
set(gcf,'Position',ss_four2three);
h_a2b3=semilogy(x_a2b3_vdx,[a2b3_VDxCox.p],'r.-');
hold on;
h_a2bInf=semilogy(x_a2bInf_vdx,[a2bInf_VDxCox.p],'k.-');
h_a2b15=semilogy(x_a2b15_vdx,[a2b15_VDxCox.p],'.-');
h_a2b10=semilogy(x_a2b10_vdx,[a2b10_VDxCox.p],'c.-');
h_a2b200=semilogy(x_a2b200_vdx,[a2b200_VDxCox.p],'m.-');
h_a2b75=semilogy(x_a2b75_vdx,[a2b75_VDxCox.p],'g.-');

semilogy(x_a2b3_vdx,repmat(0.05,size(x_a2b3_vdx)),'r--');
h_min_pval = semilogy([x_a2b3_vdx(ploc) x_a2b3_vdx(ploc)],ylim,'g--','LineWidth',2);
hold off;
xlim([0 150]);
set(gca,'fontsize',14);
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
pval_str = {['Min p-val = ', num2str(min_a2b3_VDxCox_pvals,3),' at ',...
    num2str(x_a2b3_vdx(ploc),4),' cc']};
h_lgnd=legend([h_a2b3 h_a2b10 h_a2b15 h_a2b75 h_a2b200 h_a2bInf ],...
    '$\frac{\alpha}{\beta}=3~$Gy',...
    '$\frac{\alpha}{\beta}=10~$Gy',...
    '$\frac{\alpha}{\beta}=15~$Gy',...
    '$\frac{\alpha}{\beta}=75~$Gy',...
    '$\frac{\alpha}{\beta}=200~$Gy',...
    '$\frac{\alpha}{\beta}=~$Inf',...
    'Location','Best');
set(h_lgnd,'fontsize',18);
set(h_lgnd,'interpreter','latex');

xlabel('(V_{D}) Dose [Gy]','fontsize',18); ylabel('p-value','fontsize',18);

end