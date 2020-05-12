function compute_optimal_G(scPath,scName,fcPath,startG,endG,incG,fic,ffi,simTime,dt,noiseAmp,saveMatPath)
    fc_emp = h5read(fcPath,'/CC');
    nAreas = size(fc_emp,1);
    G = startG:incG:endG;
    nIters = 1;
    fcCorrs = zeros(length(G),nIters);
    fcs = cell(length(G),nIters);
    maxfrNs = cell(length(G),nIters);
    wEI = cell(length(G),nIters);
    avgfcCorrs = zeros(length(G),1);
    for i = 1:length(G)
        Gcurr = G(i);
        avgFC = zeros(nAreas,nAreas);
        for j = 1:nIters
            [fcs{i,j}, maxfrNs{i,j}, curr, wEI{i,j}] = DMF_lesioned(scPath,scName,fic,ffi,simTime,dt,Gcurr,noiseAmp,[],[],false);
            fcCorrs(i,j) = find_corr(fcs{i,j},fc_emp);
        end
        clearvars curr;
        for k = 1:nIters
           avgFC = avgFC + fcs{i,k}; 
        end
        avgFC = avgFC/nIters;
%         [fcs{i}, frN] = DMF_lesioned(scPath,scName,fic,ffi,simTime,dt,Gcurr,noiseAmp,[],[],false);
        avgfcCorrs(i) = find_corr(avgFC,fc_emp);
        disp('saving data...');
        save(saveMatPath,'fcs','maxfrNs','avgfcCorrs','fcCorrs','wEI','nIters','G');
        disp('Done');
    end    
%     figure('Visible','off');
%     plot(G,fcCorrs,'x');
%     print('-djpeg',saveFigPath);
%     close all;
%     save(saveMatPath,'fcs','maxfrNs','fcCorrs');
end
