function generate_fc(scFnames,fcSavePath)
    for i = 1:size(scFnames,2)
       wtsFname = strsplit(scFnames{1,i},'/');
       wtsFname = wtsFname(end);
       fcFname = strcat('fc_',wtsFname);
       fc = DMF_model_bold_balance(scFnames{1,i});
       h5create(char(fullfile(fcSavePath,fcFname)),'/cc',size(fc));
       h5write(char(fullfile(fcSavePath,fcFname)),'/cc',fc);
       disp(strcat('Functional connectivity saved to ',fcSavePath,'/',fcFname));
    end
end