function DVHFromText_batch
% DVH is saved in 3 cell columns, the first stores the patient name,
% the second the original dose, original differential volume, cumulative volume from oriignal diffirential volume
% the third the resampled dose, differential volume, and cumulative volume

tic;
% parameters
pt_added=0;

dose_step=0.50;

% path and files
%pathname='Z:/elw/MATLAB/original_data/CW/CWPAIN_MASTER_DATASET_01_10_13';
pathname='Z:/elw/MATLAB/original_data/CW/CWPAIN_DATASET_01_17_13 (CW DISTANCE)';



if isunix
    %         try
    %             [xlsnum,xlstxt,xlsraw]=xlsread(pathname);
    %             save(pathname,'xlsraw');
    %         catch
    pathname=strrep(pathname,'G:','/media/SKI_G');
    load(pathname);
    %         end
else
    [~,~,xlsraw]=xlsread(pathname);
    %convert all mrns to strings
    save(pathname,'xlsraw');
end

%tmp
pathname='Z:/elw/MATLAB/original_data/CW/DVHs/';
%pathname='Z:/elw/MATLAB/original_data/CW/';
fnSuffix='_2cmExp.TXT';

% load DVH names in d
d=dir(pathname);
if isempty(d)
    disp(['!! Bad DVH location name', pathname]);
    return
end
   
f=false(size(d,1),1);
for k=1:size(d,1)
    if isempty(findstr(d(k).name,'.TXT')) || isdir(d(k).name)
        f(k)=true;
    end
end
d(f)=[];
dvh_file_names = {d.name}';
% patient last name
f=cellfun(@(x) strcmp(x,'Patient Last Name'),xlsraw(1,:)); f=find(f); f=f(1); % found the column of "Patient Last Name"
ptnames=xlsraw(2:end,f); % patient last name

% Treatment site (use to find correct DVH file for duplicate entries)
g=cellfun(@(x) strcmp(x,'Treatment Site'),xlsraw(1,:)); g=find(g); g=g(1); % found the column of "Patient Last Name"
txsites=xlsraw(2:end,g);

g=cellfun(@(x) strcmp(x,'MRN'),xlsraw(1,:)); g=find(g); g=g(1); % found the column of "MRN"
mrns=xlsraw(2:end,g);

% convert all mrn data to strings
containsNumbers = cellfun(@isnumeric,mrns);
mrns(containsNumbers) = cellfun(@num2str,mrns(containsNumbers),'UniformOutput',false);
orig_mrns=mrns;
mrns=regexprep(mrns,'-.','');
% strip

% read data from each .TXT file
DVH=cell(size(ptnames,1),4); flg_nofile=false(size(ptnames,1),1);
disp(' ');
for k=1:size(ptnames,1)
    % find the patient in the xls spreadsheet
    %f=cellfun(@(x) strcmp(x,ptnames{k}), ptnames); f=find(f);
    f=cellfun(@(x) strcmp(x,mrns{k}), mrns); f=find(f);
    
    ptname=ptnames{k};
    mrn=mrns(k);
    % remove any '-' from mrn
    dash_ind = strfind(mrn,'-');
    if ~isempty(dash_ind{1}),
        mrn = mrn{1}(1:dash_ind{1}-1);
    else
        mrn=mrn{1};
    end
    
    dvh_file_name = strcat(ptname,'_',mrn);
    
    if length(f)>1
        %disp(['repeated last names: ',ptnames{k}]);
        if length(f)>2
            disp([ptnames{k},' repeated ',num2str(length(f)),' times, quitting...']);
        end
        
        loc_str = CWLocationString(ptname,txsites{k});
        
        dvh_file_name = strcat(dvh_file_name,'_',loc_str);
    end
    % No duplicates, look up DVH with last name and MRN, strip any extra '-' first
    
    
    cur_dvh = strfind(dvh_file_names,dvh_file_name);

    dvh_name=dvh_file_names{~cellfun(@isempty,cur_dvh)};
    
    
    dvh_loc = strcat(pathname,dvh_name);
    
    disp(['Loading ',dvh_loc,'...']);
    
    % original DVH
    try
        %df=textread(strcat(pathname,ptnames{k},fnSuffix),'%s'); % read the file into cells df (data from file)
        df=textread(dvh_loc,'%s'); % read the file into cells df (data from file)
    catch ME
        disp(strcat(pathname,ptnames{k},fnSuffix)); disp(ME.message);
        %             for errcount=1:length(ME.stack)
        %                 disp(ME.stack(errcount));
        %             end
        flg_nofile(k)=true; continue;
    end
    
    DVH{k,1}=ptnames{k};
    DVH{k,4}=orig_mrns{k};
    
    df(1)=[]; % the first cell is the description text of the data, remove it
    data=cellfun(@(x) str2num(x), df);
    
    % original DVH as step function
    dose_org=data(1:4:end)/100; vol_org=data(2:4:end); % /100 to change cGy to Gy
    % cumulative dvh
    vol_cumulative=vol_org;
    %         vol_cumulative(1:end-1)=vol_cumulative(1:end-1).*diff(dose_org);
    vol_cumulative=flipud(vol_cumulative); vol_cumulative=cumsum(vol_cumulative); vol_cumulative=flipud(vol_cumulative);
    DVH{k,2}=[dose_org,vol_org,vol_cumulative];
    
    % regenerate the differentail DVH by adding new bins (which is the expected bins)
    dose_interpolation=(0:dose_step:max(DVH{k,2}(:,1))+dose_step)'; % locations of new bins
    vol_interpolation=zeros(size(dose_interpolation,1),1);
    f=ismember(dose_interpolation, dose_org); % the doses alreadyin the original dose vector should be removed in new bins
    dose_interpolation(f)=[]; vol_interpolation(f)=[];
    dose_new=[dose_org; dose_interpolation]; % add the new bins and sort it so new bins and old ones are arranged properly
    vol_new=[vol_org; vol_interpolation];
    [dose_new,g]=sort(dose_new); vol_new=vol_new(g);
    for n=1:size(dose_interpolation,1)-1 % the end point can not be interpolated, and should be zero since the last dose "bin" should correspond no volume
        f=find(dose_interpolation(n)>dose_new); % the location of interpolated dose in the DVH with combined bins, which is applied to determine the left end point (d(i-1)), d(i-1)=f(end)
        g=find(dose_interpolation(n)<dose_org); % the location of interpolated dose in the DVH with original bins, which is applied to determine the right end point (d(i)), d(i)=g(1)
        vol_new(f(end)+1) = (dose_org(g(1))-dose_interpolation(n)) / (dose_org(g(1))-dose_new(f(end))) * vol_new(f(end)); % v(j)=(d(i)-d(j))/(d(i)-d(i-1))*v(i-1)
        vol_new(f(end)) = (dose_interpolation(n)-dose_new(f(end))) / (dose_org(g(1))-dose_new(f(end))) * vol_new(f(end)); % v(i-1)=(d(j)-d(i-1))/(d(i)-d(i-1))*v(i-1)
    end
    % regenerate the differential DVH using the new bin
    dose_interpolation=(0:dose_step:max(DVH{k,2}(:,1))+dose_step)'; % locations of new bins
    vol_interpolation=zeros(size(dose_interpolation,1),1);
    for n=1:size(dose_interpolation,1)-1
        f = ( dose_interpolation(n)<=dose_new & dose_interpolation(n+1)>dose_new); % the range of new bins in the interpolated bin
        vol_interpolation(n) = sum(vol_new(f));
    end
    % cumulative dvh
    vol_cumulative=vol_interpolation;% vol_cumulative(1:end-1)=vol_cumulative(1:end-1).*diff(dose_interpolation);
    vol_cumulative=flipud(vol_cumulative); vol_cumulative=cumsum(vol_cumulative); vol_cumulative=flipud(vol_cumulative);
    DVH{k,3}=[dose_interpolation,vol_interpolation,vol_cumulative];
    
    
    %         % look up volume at each new step from original step function
    %         dose_interpolation=(0:dose_step:max(DVH{k,2}(:,1))+dose_step)'; % ((0:(ceil(max(DVH{k,1}(:,1))/dose_step)))*dose_step)';
    %         vol_interpolation=zeros(size(dose_interpolation,1),1);
    %         for n=1:size(dose_interpolation,1)
    %             f=find(dose_interpolation(n)>=dose_org);
    %             vol_interpolation(n)=vol_org(f(end));
    %         end
    %         % combine the new step function with the original one
    %         dose_new=[dose_org;dose_interpolation]; vol_new=[vol_org; vol_interpolation];
    %         [dose_new,f]=sort(dose_new); vol_new=vol_new(f);
    %         f=diff(dose_new); dose_new(f==0)=[]; vol_new(f==0)=[]; % remove duplicate records in the new step function
    %         % interpolation according to equivalent area
    %         vol_new=vol_new.*([diff(dose_new);0]); % area of each bin (step forward)
    %         vol_new(2:end)=vol_new(1:end-1); vol_new(1)=0; % area of each bin (step backward)
    %         vol_new=cumsum(vol_new); % cumulative area of the bins
    %         for n=1:size(dose_interpolation,1)
    %             f = find( dose_interpolation(n)>=dose_new ); % the cumulative interval of the new bin
    %             vol_interpolation(n)=vol_new(f(end)); % the cumulative bin
    %         end
    %         vol_interpolation(1:end-1)=diff(vol_interpolation)./diff(dose_interpolation); vol_interpolation(end)=0; % the new bin (step forward due to the diff function)
    %         % cumulative dvh
    %         vol_cumulative=vol_interpolation; vol_cumulative(1:end-1)=vol_cumulative(1:end-1).*diff(dose_interpolation);
    %         vol_cumulative=flipud(vol_cumulative); vol_cumulative=cumsum(vol_cumulative); vol_cumulative=flipud(vol_cumulative);
    %         DVH{k,3}=[dose_interpolation,vol_interpolation,vol_cumulative];
    
    % plot the result
    figure(1); a=data(1:2:end)/100; b=data(2:2:end); plot(a,b); hold on;
    stairs(DVH{k,2}(:,1),DVH{k,2}(:,2),'r');
    stairs(DVH{k,3}(:,1),DVH{k,3}(:,2),'k');
    %         a=DVH{k,2}(:,1); b=DVH{k,2}(:,2); a=[a,a]';a=a(:); a(1)=[]; b=[b,b]'; b=b(:);b(end)=[]; plot(a,b,'r');
    %         a=DVH{k,3}(:,1); b=DVH{k,3}(:,2); a=[a,a]';a=a(:); a(1)=[]; b=[b,b]'; b=b(:);b(end)=[]; plot(a,b,'k');
    %         a=dose_interpolation; b=vol_interpolation; a=[a,a]';a=a(:); a(1)=[]; b=[b,b]'; b=b(:);b(end)=[]; plot(a,b,'r');
    %         plot(dose_interpolation,vol_interpolation,'k');
    grid on; title([ptnames{k}, ' ',num2str(k)]); hold off;
    
    figure(2); stairs(DVH{k,2}(:,1),DVH{k,2}(:,3)); hold on;
    stairs(DVH{k,3}(:,1),DVH{k,3}(:,3),'r'); hold off;
    grid on; pause(0);
    pt_added=pt_added+1;
end
%     DVH(flg_nofile,:)=[]; % remove no file patients
% save result
    
    fn='Z:\elw\MATLAB\cw_analy\meta_data\CW_MASTER_DVHs.mat';
    

if isunix
    fn=strrep(fn,'G:','/media/SKI_G');
end

if 0
    disp(['Saving ',fn]);
    save(fn);
end

disp(' ');

disp(['total number of patients in .xls file: ',num2str(pt_added)]);
disp(['total number of patients with the DVH files: ', num2str(sum(~flg_nofile))]);
disp(['total number of patients without the DVH files: ', num2str(sum(flg_nofile))]);
toc;