function DataFromCW % extracting for MD analysis
tic;
% parameters
    dvhdef={'2cmExp';'Colorado'};
    CG = cell(length(dvhdef),1);

% load data
    for k = 1:length(dvhdef)
        fn=['G:/MSKCC/Andy/Ken/tom/CW_',dvhdef{k},'.mat'];
        if isunix
            fn=strrep(fn,'G:','/media/SKI_G');
        end
        load(fn,'CGobj');
        CG{k} = CGobj;
    end

% extract V20 and V30
    Dose = [20; 30];
    for k = 1:length(dvhdef)
        % extract data
        vx = zeros(CGobj.mNumInGrp,length(Dose));
        for n = 1:length(Dose)
            for m = 1:CGobj.mNumInGrp
                vx(m,n) = CG{k}.mGrp(m).fVolAtDose(Dose(n));
            end
        end
        % write to excel file
        d = cell(CGobj.mNumInGrp+1,length(Dose)+1);
        d{1,1} = 'Pt Code';
        for n = 1:length(Dose)
            d{1,n+1} = ['V',num2str(Dose(n))];
        end
        d(2:end,1) = {CGobj.mGrp.mID};
        d(2:end,2:length(Dose)+1) = num2cell(vx);

        fn=['G:/MSKCC/Andy/Ken/tom/CW_Vx_',dvhdef{k},'.xls'];
        if isunix
            fn=strrep(fn,'G:','/media/SKI_G');
        end
        xlswrite(fn,d,dvhdef{k});
    end
toc;
end