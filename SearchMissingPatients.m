function SearchMissingPatients
tic;
% parameters;
    flg_2cmExp=true;
% path and files
    pathname='G:/MSKCC/Andy/Ken/meta/3_28_2010_CW_final_cohort';
    if isunix
        pathname=strrep(pathname,'G:','/media/SKI_G');
    end
    try
        load([pathname,'.mat']);
    catch
        warning('off','MATLAB:xlsreadold:Truncation');
        [~,~,xlsraw]=xlsread(strcat(pathname,'.xls'));
        warning('on','MATLAB:xlsreadold:Truncation');
        save(pathname,'xlsraw');
    end
    
    if flg_2cmExp
        pathname='//Pensmph6/MpcsResearch1/JacksonA/CWDVHs/HFX_Chestwall/2cmExpansion/';
        fnSuffix='_2cmExp.TXT';
    else
        pathname='//Pensmph6/MpcsResearch1/JacksonA/CWDVHs/HFX_Chestwall/Colorado/';
        fnSuffix='_Colorado.TXT';
    end
    if isunix
        pathname=strrep(pathname,'H:','/media/SKI_H');
    end

% files in the path
    d=dir(pathname);
    f=false(size(d,1),1);
    for k=1:size(d,1)
        if isempty(findstr(d(k).name,'.TXT')) || isdir(d(k).name)
            f(k)=true;
        end
    end
    d(f)=[];

% patient last name in the .xls file
    f=cellfun(@(x) strcmp(x,'Patient Last Name'),xlsraw(1,:)); f=find(f); f=f(1); % found the column of "Patient Last Name"
    ptnames=xlsraw(2:end,f); % patient last name

% check patient names in the spread sheet using the file names
    disp(' ');
    for k=1:size(d,1)
        pt=strrep(d(k).name,fnSuffix,'');
        f=cellfun(@(x) strcmp(x,pt), ptnames);
        if all(~f)
            disp(d(k).name);
        end
    end
toc;