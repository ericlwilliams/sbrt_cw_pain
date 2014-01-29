function ChestWallPainAlpha2BetaDv
tic;

fig_loc = 'Z:/elw/MATLAB/cw_analy/figures/latest/';

fp = 'Z:\elw\MATLAB\cw_analy\meta_data\';
%fp = 'C:\Documents and Settings\williae1\cw_meta_data\';

fig_num=1;
if isunix
    fp=strrep(fp,'G:','/media/SKI_G');
end
%fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_ratio0.mat'};
fn = {'MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};
%fn = {'RIMNER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat'};

CGobj = cell(length(fn),1);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

% load data

load(strcat(fp,fn{1}),'CGobj_current');
CGobj = CGobj_current;
CGgrp = CGobj.mGrp;
n_grp = CGobj.mNumInGrp;


%max_vol = max(max(CGgrp.mVolCum));

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

a2b_dx = zeros(length(a2b),1);
a2b_logl = -inf(length(a2b),1);

fig=figure(1);
set(gcf,'Position',ss_four2three);
hold on;

colors = varycolor(length(a2b));
% colors = cell(30,1);
% cur_color = [1:-0.05:0.5]';
% for l=1:length(cur_color);
%     colors{l}=[cur_color(l) 0 0];
%     colors{l+10}=[0 cur_color(l) 0];
%     colors{l+20}=[0 0 cur_color(l)];
% end

dx_logls = cell(length(a2b),1);
med_dx = -inf(length(a2b),1);

for i=1:length(a2b)
    
    cur_a2b = a2b(i);
    
    % loop over V1 - V100
    max_dx=200; %200 cc
    
    cur_dx_logls = -inf(max_dx,1);
    cur_med_dx = -inf(max_dx,1);
    
    
    for k=1:max_dx
        volume = k;
        
        dv=zeros(n_grp,1); % dose at volume
        dv(:)=0;
    
                %loop of patients, calculate new dosebins and vXs
        for j=1:n_grp
            
            num_fx = CGgrp(j).mFxNum;
            cur_dosebins = CGgrp(j).mDoseBins_org;
            vols = CGgrp(j).mVolCum;
            if cur_a2b<Inf
                new_dosebins = cur_dosebins.* ...
                    ((cur_a2b + (cur_dosebins/num_fx))/...
                    (cur_a2b + 2));
            else
                new_dosebins = cur_dosebins;
            end
                        
            %dv(j) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(fdose) );
            
            %f = find( new_dosebins <= dose ); f=f(end); % f won't be empty by definition
            f = find( vols > volume );
            if isempty(f)
                cv(j) = 0;
                continue;
            end
                f=f(end); % f won't be empty by definition
            if vols(f) < vols(1)
                %dv(j) = interp1( [new_dosebins(f); new_dosebins(f+1)], [vols(f); vols(f+1)], dose );
                dv(j) = interp1( [vols(f);vols(f+1)], [new_dosebins(f); new_dosebins(f+1)], volume);
            else
                dv(j) = 0;
            end
            
     
        end

        if sum(dv)==0
            break;
        end
        %[~,logl,~,~]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
        [~,logl,~,~]=coxphfit(dv,compdate,'censoring',flgcensor);
        cur_dx_logls(k) = logl;
        cur_med_dx(k) = median(dv);
    end


    
     p=plot(1:max_dx,cur_dx_logls);
    set(p,'Color',colors(i,:))
    [a2b_logl(i),a2b_dx(i)] = max(cur_dx_logls);
    
    med_dx(i) = cur_med_dx(a2b_dx(i));
    dx_logls{i}=cur_dx_logls;
end

 
%print_fig(fig,fig_loc,'a2b_dx_fits','png');

fn=['ChestWall_a2b_dv.mat'];
            
save(strcat(fp,fn),'a2b','a2b_logl','a2b_dx','dx_logls',...
    'compdate','flgcensor');
end