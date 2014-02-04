function ChestWallPainAnalysis
tic;

save_result = true;

% parameters
dvhdef={'DVHs'};
fxnum={-1};
%fxnum={3 4 5}; % -1 all fractions, [n1 n2 ...] analyze patients with fractions of n1, n2, ... treatment planning only
%beta2alpha=[1/2.1];
%beta2alpha=[0];
beta2alpha=[Inf];
a2b_corr = 'NTD';
dosestep=1;
volstep=1; % volume step in cc
timestep=3; % time step in month

% patient info.
% load patient info stored in the spread sheet
fn='Z:\elw\MATLAB\original_data\CW\CWPAIN_DATASET_01_17_13 (CW DISTANCE)';


if isunix
    fn=strrep(fn,'G:','/media/SKI_G');
end
try
    load(fn,'xlsraw');
catch
    error(['Error loading CW PAIN DATASET']);
    fn1 = strrep(fn,'tom','meta');
    [~,~,xlsraw]=xlsread([fn1,'.xlsx'],'sheet1');
    save(fn,'xlsraw');
end

% pick up related data
PtInfo = classDataFromXls();
PtInfo.mXlsRaw = xlsraw;
% MRN
%PtInfo.mLabel = 'Patient Last Name';
PtInfo.mLabel = 'MRN';
PtInfo = PtInfo.fExtractColData();
flgptcode = PtInfo.mFlg;
ptcode = PtInfo.mData;
% gender
PtInfo.mLabel = 'Sex';
PtInfo = PtInfo.fExtractColData();
flggender = PtInfo.mFlg;
ptgender = PtInfo.mData;
ptgender = strcmp(ptgender,'Male');

% Number of Fractions
PtInfo.mLabel = 'Number of Fractions';
PtInfo = PtInfo.fExtractColData();
flgfx = PtInfo.mFlg;
fx = zeros(size(flgfx));
fx(flgfx) = cell2mat(PtInfo.mData(flgfx));

% KPS
PtInfo.mLabel = 'KPS at dx';
PtInfo = PtInfo.fExtractColData();
flgkps = PtInfo.mFlg;
kps = zeros(size(flgkps));
kps(flgkps) = cell2mat(PtInfo.mData(flgkps));

% Distance to CW
PtInfo.mLabel = 'Distance to CW (cm)';
PtInfo = PtInfo.fExtractColData();
flgcm2cw = PtInfo.mFlg;
cm2cw = zeros(size(flgcm2cw));
cm2cw(flgcm2cw) = cell2mat(PtInfo.mData(flgcm2cw));


% Delivered dose
PtInfo.mLabel = 'Total Dose (cGy)';
PtInfo = PtInfo.fExtractColData();
flgtx = PtInfo.mFlg;
tx = zeros(size(flgtx));
tx(flgtx) = cell2mat(PtInfo.mData(flgtx));
% Date of Birth
PtInfo.mLabel = 'Date of Birth';
PtInfo = PtInfo.fExtractColData();
flgdatebirth = PtInfo.mFlg;
datebirth = zeros(size(flgdatebirth));
f = PtInfo.mData(flgdatebirth);
datebirth(flgdatebirth) = datenum(f);
% Date of Birth


% IGRT End Date
PtInfo.mLabel = 'IGRT End Date';
PtInfo = PtInfo.fExtractColData();
flgdateIGRT = PtInfo.mFlg;
dateIGRT = zeros(size(flgdateIGRT));
f = PtInfo.mData(flgdateIGRT);
dateIGRT(flgdateIGRT) = datenum(f);
% IGRT Start Date
PtInfo.mLabel = 'IGRT Start Date';
PtInfo = PtInfo.fExtractColData();
flgdateStartIGRT = PtInfo.mFlg;
dateStartIGRT = zeros(size(flgdateStartIGRT));
f = PtInfo.mData(flgdateStartIGRT);
dateStartIGRT(flgdateStartIGRT) = datenum(f);

% age (at start of treatment)
flgAge = (flgdateStartIGRT & flgdatebirth);
age = (dateStartIGRT-datebirth)./365;

% baseline bmi
PtInfo.mLabel = 'Baseline BMI';
PtInfo = PtInfo.fExtractColData();
flgbmi = PtInfo.mFlg;
data = PtInfo.mData(flgbmi);
idx = cellfun(@ischar,data);
data(idx)={[-Inf]};% patients without bmi

bmi = [0;cell2mat(data)];




PtInfo.mLabel = 'Date of last follow-up';
PtInfo = PtInfo.fExtractColData();
flgdatefu = PtInfo.mFlg;
datefu = zeros(size(flgdatefu));
datefu(flgdatefu) = datenum(PtInfo.mData(flgdatefu));

% complication grade 2
PtInfo.mLabel = 'CW Pain date of Grade 2 (ROB MUTTER CLASS)';
%PtInfo.mLabel = 'CW Pain date of Grade 2 (RECENT CLASS)';% Rimner
%PtInfo.mLabel = 'Rib fracture date';

PtInfo = PtInfo.fExtractColData();
datepain = inf(size(PtInfo.mData));
flgcensor = true(size(datefu));

f = PtInfo.mFlg;
datepain(f) = datenum(PtInfo.mData(f));
flgcensor(f) = false;



%%
% relapse date
%PtInfo.mLabel = 'Local Failure';
PtInfo.mLabel = 'CW Pain resolved? (0=No, 1=Yes)';
PtInfo = PtInfo.fExtractColData();
f = PtInfo.mFlg;
flgfailure = false(size(f));
flgfailure(f) = cell2mat(PtInfo.mData(f));
%PtInfo.mLabel = 'Local Failure Date';
%% Should be 'CW Pain recurrence date'!
%PtInfo.mLabel = 'Date of last follow-up';
PtInfo.mLabel = 'CW Pain recurrence date';
PtInfo = PtInfo.fExtractColData();
flgrelapse = PtInfo.mFlg;
%datefailure = inf(size(flgrelapse));
%datefailure(flgrelapse&flgfailure) = datenum(PtInfo.mData(flgrelapse&flgfailure));
%flgrelapse = flgrelapse&f;


% death date
%PtInfo.mLabel = 'Surv alive 0, dead 1';
PtInfo.mLabel = 'Survival (0 = Alive, 1 = Dead)';
PtInfo = PtInfo.fExtractColData();
f1 = PtInfo.mFlg;
flgdeath = false(size(f1));
flgdeath(f1) = cell2mat(PtInfo.mData(f1));
%PtInfo.mLabel = 'death Date';
PtInfo.mLabel = 'Survival Date/Date of Death';
PtInfo = PtInfo.fExtractColData();
f = PtInfo.mFlg;
datedeath = inf(size(f));
datedeath(f&flgdeath) = datenum(PtInfo.mData(f&flgdeath));
flgdeath = f1&f;



% common part of those data

%%TMP: not checking/loading recurance date
flg = flgptcode & flgfx & flgkps & flgcm2cw & flgtx & flgdateIGRT &...
    flgdatefu & flgdatebirth & flggender &...
    flgdeath & flgAge & flgbmi & flgdateStartIGRT;

%flg = flgptcode & flgfx & flgtx & flgdatebirth & flgdateIGRT & flgdatefu & flgrelapse & flggender & flgdeath;


ptcode=ptcode(flg);
% convert all mrn data to strings
containsNumbers = cellfun(@isnumeric,ptcode);
ptcode(containsNumbers) = cellfun(@num2str,ptcode(containsNumbers),...
    'UniformOutput',false);
fx=fx(flg);
kps=kps(flg);
cm2cw=cm2cw(flg);
tx=tx(flg);
datebirth=datebirth(flg);
dateIGRT=dateIGRT(flg);
dateStartIGRT=dateStartIGRT(flg);
datefu=datefu(flg);
datepain=datepain(flg);
flgcensor=flgcensor(flg);
%datefailure = datefailure(flg);
ptgender = ptgender(flg);
age = age(flg);
bmi = bmi(flg);
datedeath = datedeath(flg);

% patient DVH and objects
for m=1:length(dvhdef) % iterate with each definition
    % load dvh info
    fn=['Z:\elw\MATLAB\cw_analy\meta_data\CW_MASTER_',dvhdef{m}];
    
    if isunix
        fn=strrep(fn,'G:','/media/SKI_G');
    end
    load(fn,'DVH');
    
    % match patients DVH and .xls
    [f1,g1]=ismember(ptcode,DVH(:,4));
    [f2,g2]=ismember(DVH(:,4),ptcode);
    if length(unique(g2(g2~=0)))~=length(find(g2))% || length(unique(g2~=0))~=length(find(g2))
        error('The patient ids are not unique');
    end
    if ~( all(f1) && all(f2) )
        disp('The patients in the spread sheet do not match with those in DVH curves, take the common part');
        % use common part only
        DVH=DVH(f2,:);
        ptcode=ptcode(f1); fx=fx(f1); tx=tx(f1); kps=kps(f1);
        datebirth=datebirth(f1); cm2cw=cm2cw(f1);
        dateStartIGRT=dateStartIGRT(f1);dateIGRT=dateIGRT(f1);
        datefu=datefu(f1); datepain=datepain(f1); flgcensor=flgcensor(f1);
        [~,g1]=ismember(ptcode,DVH(:,1));
    end
    
    % combine info from DVH and .xls to form complicationgroup object
    CIobjs = classOutcomeIndividual();
    CIobjs.mBeta2AlphaCorrection = a2b_corr;
    %CIobjs.mLymanN = 10.^(-1:0.1:1)';
    CIobjs = repmat(CIobjs,[size(DVH,1),1]);
    for n = 1:size(DVH,1)
        CIobjs(n,1).mID=ptcode{n};
        CIobjs(n,1).mGender = ptgender(n);
        CIobjs(n,1).mAgeAtTx = age(n);
        CIobjs(n,1).mBMI = bmi(n);
        CIobjs(n,1).mFxNum=fx(n);
        CIobjs(n,1).mKPS=kps(n);
        CIobjs(n,1).mDistanceToChestWall=cm2cw(n);
        CIobjs(n,1).mDoseTx=tx(n);
        CIobjs(n,1).mDosePerFx=tx(n)/fx(n);
        CIobjs(n,1).mDateBirth = datebirth(n);
        CIobjs(n,1).mDateDeath = datedeath(n);
        CIobjs(n,1).mDateBaseline = dateIGRT(n);
        CIobjs(n,1).mDateStartTx = dateStartIGRT(n);
        CIobjs(n,1).mDateComp = datepain(n);
        CIobjs(n,1).mDateLastFollowup = datefu(n);
        %CIobjs(n,1).mDateRelapse = datefailure(n);
        CIobjs(n,1).mFlgCensor = flgcensor(n);
        CIobjs(n,1).mDoseBins_org=DVH{g1(n),2}(:,1);
        CIobjs(n,1).mVolDiff=DVH{g1(n),2}(:,2);
        CIobjs(n,1).mVolCum=DVH{g1(n),2}(:,3);
        
    end
    
    if any([CIobjs.mFxNum]==0)
        error('Fraction number is zeros');
    end
    CGobj_org = classOutcomeAnalysis();
    %CGobj_org.mLymanN = 10.^(-1:0.1:1)';
    CGobj_org.mStepDose = dosestep;
    CGobj_org.mStepVol = volstep;
    CGobj_org.mStepTime = timestep;
    

    
    CGobj_org = CGobj_org.fAddPatient(CIobjs);
%    CGobj_org = CGobj_org.fCalculateEUD();
    
    
    % per fraction category
    for n=1:length(fxnum)
        % pick up patient with the wanted fraction numbers
        if any(fxnum{n}==-1)
            f1=true(CGobj_org.mNumInGrp,1);
        else % retrieve patients with exact fractions only
            f1=arrayfun(@(x) any(x==fxnum{n}),[CGobj_org.mGrp.mFxNum]);
        end
        CGobj_current = CGobj_org;
        CGobj_current.mGrp = CGobj_current.mGrp(f1);
        CGobj_current.mNumInGrp = size(CGobj_current.mGrp,1);
        
        % up from 2 to prevent hanging, could specify only with 'new'
        %CGobj_current.mCoxMinSize = 20;
        
        % for each alpha/beta compute a series of complication table
        disp(['Running KM and CPH analyses...']);
        for u = 1:length(beta2alpha)
            % adjust beta to alpha ratio

            CGobj_current.mBeta2Alpha = beta2alpha(u);
            %CGobj_current = CGobj_current.fLinearQuartraticCorrection();
             
             
        
            disp(['LogRankCM2CW_DVH...']); % V_{x} with empty patients set to zero
             CGobj_current = CGobj_current.fLogRankCM2CW_DVH();
            
             disp(['CoxModel_DVH...']);
            CGobj_current = CGobj_current.fCoxModel_DVH();
            
              % calculate log-rank pval
             disp(['LogRankV30BMITx_DVH...']); % V_{x} with empty patients set to zero
             CGobj_current = CGobj_current.fLogRankV30BMITx_DVH();
%             
             disp(['LogRankV39BMI_DVH...']); % V_{x} with empty patients set to zero
             CGobj_current = CGobj_current.fLogRankV39BMI_DVH();
%         
             disp(['LogRankV30BMI_DVH...']); % V_{x} with empty patients set to zero
             CGobj_current = CGobj_current.fLogRankV30BMI_DVH();
%             
             disp(['LogRankBMI_DVH...']); % V_{x} with empty patients set to zero
             CGobj_current = CGobj_current.fLogRankBMI_DVH();
             
            disp(['LogRankVDx_DVH...']); % V_{x} with empty patients set to zero
            CGobj_current = CGobj_current.fLogRankVDx_DVH();
            
            % analysis 
            
            disp(['OveralSurvivalCurve...']);
            CGobj_current = CGobj_current.fOverallSurviveCurve();
            disp(['OveralCompCurve...']);
            CGobj_current = CGobj_current.fOverallCompCurve();

           
         
            
%             disp(['LogRankDx_DVH...']); % D_{x} with empty patients excluded
%             CGobj_current = CGobj_current.fLogRankDx_DVH();
%             

          
         
            % contains same info as above 
            %disp(['LogRankDVx_DVH...']);% D_{x} with empty patients set to zero
            %CGobj_current = CGobj_current.fLogRankDVx_DVH();
            
            % save result
            CGstrct = ObjToStruct(CGobj_current);
            %fn=['Z:\elw\MATLAB\cw_analy\meta_data\MUTTER_ChestWall_Cox_DiVj_',dvhdef{m},'_fx',num2str(fxnum{n}(1)),'_ratio',num2str(beta2alpha(u)),'.mat'];
            %fn=['Z:\elw\MATLAB\cw_analy\meta_data\MUTTER_EARLY_ChestWall_Cox_DiVj_',dvhdef{m},'_fx',num2str(fxnum{n}(1)),'_ratio',num2str(beta2alpha(u)),'.mat'];
            %fn=['Z:\elw\MATLAB\cw_analy\meta_data\RIMNER_EARLY_ChestWall_Cox_DiVj_',dvhdef{m},'_fx',num2str(fxnum{n}(1)),'_ratio',num2str(beta2alpha(u)),'.mat'];
            %fn=['Z:\elw\MATLAB\cw_analy\meta_data\RIMNER_LATE_ChestWall_Cox_DiVj_',dvhdef{m},'_fx',num2str(fxnum{n}(1)),'_ratio',num2str(beta2alpha(u)),'.mat'];
            
            fn=['C:\Documents and Settings\williae1\cw_meta_data\MUTTER_MASTER_ChestWall_Cox_DiVj_',dvhdef{m},'_fx',num2str(fxnum{n}(1)),'_a2b',num2str(1/beta2alpha(u)),'.mat'];
           
            if isunix
                fn=strrep(fn,'G:','/media/SKI_G');
            end
            
            if save_result
                save(fn,'CGobj_current','CGstrct');
                disp(fn);
            end
        end
    end
end
toc;