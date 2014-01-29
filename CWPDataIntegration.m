function CWPDataIntegration
tic;
% parameters
    dvhdef={'2cmExp';'Colorado'};

    b2a = 0;

% load patient info stored in the spread sheet
	fn='Z:/Fan/Andy/Ken/meta/20100424CWFinalCohort';
    if isunix
        fn=strrep(fn,'G:','/media/SKI_G');
    end
    try
        load(fn,'xlsraw');
    catch
        [~,~,xlsraw]=xlsread([fn,'.xlsx'],'sheet1');
        save(fn,'xlsraw');
    end

% pick up related data from the spread sheet
    PtInfo = classDataFromXls();
    PtInfo.mXlsRaw = xlsraw;
    % identity
    PtInfo.mColName = 'Patient Last Name';
    PtInfo = PtInfo.fExtractColData();
    flgptcode = PtInfo.mFlgDataRows;
    ptcode = PtInfo.mColData;
    % Number of Fractions
    PtInfo.mColName = 'Number of Fractions';
    PtInfo = PtInfo.fExtractColData();
    flgfx = PtInfo.mFlgDataRows;
    fx = zeros(size(flgfx));
    fx(flgfx) = cell2mat(PtInfo.mColData(flgfx));
    % Date of Birth
    PtInfo.mColName = 'Date of Birth';
    PtInfo = PtInfo.fExtractColData();
    flgdatebirth = PtInfo.mFlgDataRows;
    datebirth = zeros(size(flgdatebirth));
    f = PtInfo.mColData(flgdatebirth);
    datebirth(flgdatebirth) = datenum(f);
    % IGRT Start Date
    PtInfo.mColName = 'IGRT Start Date';
    PtInfo = PtInfo.fExtractColData();
    flgdateIGRT = PtInfo.mFlgDataRows;
    dateIGRT = zeros(size(flgdateIGRT));
    f = PtInfo.mColData(flgdateIGRT);
    dateIGRT(flgdateIGRT) = datenum(f);
    % Last Follow Up Date
    PtInfo.mColName = 'Surv Date';
    PtInfo = PtInfo.fExtractColData();
    flgdatefu = PtInfo.mFlgDataRows;
    datefu = zeros(size(flgdatefu));
    datefu(flgdatefu) = datenum(PtInfo.mColData(flgdatefu));
    % complication grade 3
    PtInfo.mColName = 'Pain Grade 3';
    PtInfo = PtInfo.fExtractColData();
    datepain = inf(size(PtInfo.mColData));
    flgcensor = true(size(datefu));
    f = PtInfo.mFlgDataRows;
    datepain(f) = datenum(PtInfo.mColData(f));
    flgcensor(f) = false;
    % complication grade 2
    PtInfo.mColName = 'Pain Grade 2';
    PtInfo = PtInfo.fExtractColData();
    f = PtInfo.mFlgDataRows;
    datepain(f) = datenum(PtInfo.mColData(f));
    flgcensor(f) = false;
    % relapse date
    PtInfo.mColName = 'Local Failure';
    PtInfo = PtInfo.fExtractColData();
    f = PtInfo.mFlgDataRows;
    flgfailure = false(size(f));
    flgfailure(f) = cell2mat(PtInfo.mColData(f));
    PtInfo.mColName = 'Local Failure Date';
    PtInfo = PtInfo.fExtractColData();
    flgrelapse = PtInfo.mFlgDataRows;
    datefailure = inf(size(flgrelapse));
    datefailure(flgrelapse&flgfailure) = datenum(PtInfo.mColData(flgrelapse&flgfailure));
    flgrelapse = flgrelapse&f;
    % gender
    PtInfo.mColName = 'Sex';
    PtInfo = PtInfo.fExtractColData();
    flggender = PtInfo.mFlgDataRows;
    ptgender = PtInfo.mColData;
    % death date
    PtInfo.mColName = 'Surv alive 0, dead 1';
    PtInfo = PtInfo.fExtractColData();
    f1 = PtInfo.mFlgDataRows;
    flgdeath = false(size(f1));
    flgdeath(f1) = cell2mat(PtInfo.mColData(f1));
    PtInfo.mColName = 'death Date';
    PtInfo = PtInfo.fExtractColData();
    f = PtInfo.mFlgDataRows;
    datedeath = inf(size(f));
    datedeath(f&flgdeath) = datenum(PtInfo.mColData(f&flgdeath));
    flgdeath = f1&f;
    
    % common part of those data
    flg = flgptcode & flgfx & flgdatebirth & flgdateIGRT & flgdatefu & flgrelapse & flggender & flgdeath;
    ptcode=ptcode(flg);
    fx=fx(flg);
    datebirth=datebirth(flg);
    dateIGRT=dateIGRT(flg);
    datefu=datefu(flg);
    datepain=datepain(flg);
    flgcensor=flgcensor(flg);
    datefailure = datefailure(flg);
    ptgender = ptgender(flg);
    datedeath = datedeath(flg);
        
% patient DVH and objects
    CGobj_org=classEndPointGroup();
    CGobj_org.mBeta2Alpha = b2a;
    for m=1:length(dvhdef) % iterate with each definition
        % load dvh info
        fn=['G:/MSKCC/Andy/Ken/meta/DVH_',dvhdef{m}];
        if isunix
            fn=strrep(fn,'G:','/media/SKI_G');
        end
        load(fn,'DVH');
        
        % match patients DVH and .xls
        [f1,g1]=ismember(ptcode,DVH(:,1));
        [f2,g2]=ismember(DVH(:,1),ptcode);
        if length(unique(g2(g2~=0)))~=length(find(g2))% || length(unique(g2~=0))~=length(find(g2))
            error('The patient ids are not unique');
        end
        if ~( all(f1) && all(f2) )
            disp('The patients in the spread sheet do not match with those in DVH curves, take the common part');
            % use common part only
            DVH=DVH(f2,:);
            ptcode=ptcode(f1); fx=fx(f1); datebirth=datebirth(f1); dateIGRT=dateIGRT(f1); datefu=datefu(f1); datepain=datepain(f1); flgcensor=flgcensor(f1);
            [~,g1]=ismember(ptcode,DVH(:,1));
        end
        
        % combine info from DVH and .xls to form complicationgroup object
        CIobjs = classEndPointIndividual();
        CIobjs = repmat(CIobjs,[size(DVH,1),1]);
        for n = 1:size(DVH,1)
            CIobjs(n,1).mID=ptcode{n};
            CIobjs(n,1).mDoseBins_org=DVH{g1(n),2}(:,1);
            CIobjs(n,1).mVolDiff=DVH{g1(n),2}(:,2);
            CIobjs(n,1).mVolCum=DVH{g1(n),2}(:,3);
            CIobjs(n,1).mFxNum=fx(n);
            %                 CIobjs(n,1).DoseStep = dosestep;
            %                 CIobjs(n,1).VolStep = volstep;
            CIobjs(n,1).mDateBirth = datebirth(n);
            CIobjs(n,1).mDateBaseline = dateIGRT(n);
            CIobjs(n,1).mDateComp = datepain(n);
            CIobjs(n,1).mDateLastFollowup = datefu(n);
            CIobjs(n,1).mFlgCensor = flgcensor(n);
            CIobjs(n,1).mDateRelapse = datefailure(n);
            CIobjs(n,1).mGender = ptgender(n);
            CIobjs(n,1).mDateDeath = datedeath(n);
        end
        if any([CIobjs.mFxNum]==0)
            error('Fraction number is zeros');
        end
        CGobj = CGobj_org.fAddPatient(CIobjs);
        
        % save result
        CGstrct = ObjToStruct(CGobj);
        fn=['G:/MSKCC/Andy/Ken/tom/CW_',dvhdef{m},'.mat'];
        if isunix
            fn=strrep(fn,'G:','/media/SKI_G');
        end
        save(fn,'CGobj','CGstrct');
        disp(fn);
    end
toc;
end