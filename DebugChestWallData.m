function DebugChestWallData
tic;

[~,~,comparision_data]=xlsread('Z:\elw\MATLAB\original_data\CW\CWPAIN_DATASET_OLDNEW_RP_COMPARISION.xlsx');

pt_name = comparision_data(2:end,2);

mutter_rp2 = comparision_data(:,12);
mutter_rp2 = mutter_rp2(2:end);

mutter_rp2_outcome = ones(size(mutter_rp2));
mutter_rp2_date = zeros(size(mutter_rp2));
for i=1:length(mutter_rp2_outcome);

    if isempty(mutter_rp2{i}),
        mutter_rp2_outcome(i)=0;
    else
        mutter_rp2_date(i)=datenum(mutter_rp2{i});
    end
end


rimner_rp2 = comparision_data(:,17);
rimner_rp2 = rimner_rp2(2:end);

rimner_rp2_outcome = ones(size(rimner_rp2));
rimner_rp2_date = zeros(size(rimner_rp2));
for i=1:length(rimner_rp2_outcome);

    if isnan(rimner_rp2{i}),
        rimner_rp2_outcome(i)=0;
    else
        rimner_rp2_date(i)=datenum(rimner_rp2{i});
    end
end



mutter_fu_date = comparision_data(:,19);
mutter_fu_date = mutter_fu_date(2:end);
mutter_fu_date= datenum(mutter_fu_date);

%% find problems



rim_minus_mutt = rimner_rp2_outcome-mutter_rp2_outcome;

comp_date_diff_n=0;
for j=1:length(rim_minus_mutt)
    cur_diff = rim_minus_mutt(j);
   
    %find problems
    if (cur_diff == 0),
        if (rimner_rp2_outcome(j)==1),
            %both have gd2 comp, check that dates are the same
            if rimner_rp2_date(j)~=mutter_rp2_date(j),
                comp_date_diff_n=comp_date_diff_n+1;
            
%                 disp(['-']);
%                 disp(['Comp dates differ for:']);
%                 disp([pt_name{j},...
%                     ' -> R: ',datestr(rimner_rp2_date(j)),...
%                     ' -> M: ',datestr(mutter_rp2_date(j))]);
            end
        end
    elseif (cur_diff > 0), % new data has complication, old doesn't
        %make sure that rimner_rp2_date > mutter_date_last_fu
        if rimner_rp2_date(j)<mutter_fu_date(j),
%             disp(['%%'])
%             disp(['New comp, unmatched']);
%             disp([pt_name{j},...
%                 ' R: ',datestr(rimner_rp2_date(j)),...
%                 ' M_FU: ',datestr(mutter_fu_date(j))]);
        end
    else % mutter comp no rimner compe
        disp(['$$'])
        disp(['OLD comp, unmatched:']);
        disp([pt_name{j}]);
    end
    
end
disp(['*******']);
disp([num2str(comp_date_diff_n), ' patients have different 2gd complication dates']);
toc;
end