startTime = tic;
clear;
scPaths = {'../../data/connectivity_998/Hagmann_998.hdf5'};
fcPaths = {'../../data/connectivity_998/Hagmann_998.hdf5'};
resultsPaths = {'../../results/998/simhc_hag998_resampled_G0.2_dt1_sigma0.001_8mins.hdf5'};
scName = '/C_resampled';
nAreas = 998;
nItersPerSim = 1;
simTime = 8*60*1000;
fic = true;
ffi = false;
dt = 1;
G = 0.5;
noiseAmp = 0.001;
fc_emp = h5read(fcPaths{1},'/CC');
for i = 1:length(scPaths)
    avgCorr = 0;
    maxFCCorr = 0;
    scPath = scPaths{i};
    resultsPath = resultsPaths{i};
    rmCmd = strcat(['rm ',resultsPath]);
    if(~system(strcat(['test -e ',resultsPath])))
        system(rmCmd);
    end

    for j = 1:nItersPerSim
        [fc, frN, curr, wEI]=DMF_lesioned(scPath,scName,fic,ffi,simTime,dt,G,noiseAmp,[],[],false);
        if(j == 1)
            avgFC = fc;
        else
            avgFC = avgFC + fc;
        end
        fcCorr = find_corr(fc,fc_emp);
        disp(strcat(['fcCorr: ',num2str(fcCorr)]));
        avgCorr = avgCorr + fcCorr;
        if(maxFCCorr < fcCorr)
            maxFCCorr = fcCorr;
            maxFC = fc;
            maxFCfrN = frN;
        end 
    end
    avgFC = avgFC / nItersPerSim;
    avgCorr = avgCorr/nItersPerSim;
    avgFCCorr = find_corr(avgFC,fc_emp);
    frNmean = mean(maxFCfrN);
    frNmax = max(maxFCfrN);
    h5create(resultsPath,'/CC',[nAreas nAreas]);
    h5write(resultsPath,'/CC',maxFC);
    h5create(resultsPath,'/avgCC',[nAreas nAreas]);
    h5write(resultsPath,'/avgCC',avgFC);
    h5create(resultsPath,'/wEI',size(wEI));
    h5write(resultsPath,'/wEI',wEI);
    h5create(resultsPath,'/frNmean',nAreas);
    h5write(resultsPath,'/frNmean',frNmean);
    h5create(resultsPath,'/frNmax',nAreas);
    h5write(resultsPath,'/frNmax',frNmax);


    disp(strcat(['Avg. FC correlatoin for ',scPath,': ',num2str(avgFCCorr)]));
    disp(strcat(['Avg. correlatoin for ',scPath,': ',num2str(avgCorr)]));
    disp(strcat(['Max. FC correlatoin for ',scPath,': ',num2str(maxFCCorr)]));
%     disp(strcat(['Saved FC of ',scPath,' as: ',resultsPath]));
end
toc(startTime);
