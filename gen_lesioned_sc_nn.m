function [lSC] = gen_lesioned_sc_nn(scPath,centresPath,scName,centresName,lesionIdx,lesionPercent,saveDir)
    disp('WARNING:This function assumes that the first half rows are areas in right hemisphere and last half are areas in left hemisphere');
%     hag_data = load(scPath);
    lSC = h5read(scPath,scName);
%     lSC = hag_data.CIJ_fbden_average;
    nAreas = size(lSC,1);
    centres = h5read(centresPath,centresName);
%     centres_66 = h5read('../../data/connectivity_66/Hagmann_66.hdf5','/Talairach');
%     centres_998 = h5read(scPath,'/Talairach');
%     regionNames = strsplit(fileread('../../data/connectivity_68/fs_region_abbrev_sort.txt'),'\n');
%     regionNames = h5read(scPath,'/ROI_lbls');
%     lesionIdx = lesionIdxs(1);
%     lesionName = regionNames{lesionIdx};
    nNeighbors = floor(nAreas*lesionPercent) + 1;
%     if(strcmp(hemisphere,'right'))
%         startIdx = 1;
%         endIdx = nAreas/2;
%     else
%         startIdx = nAreas/2 + 1;
%         endIdx = nAreas;
%     end 
    
    if(lesionIdx <= nAreas/2)
        startIdx = 1;
        endIdx = nAreas/2;
    else
        startIdx = nAreas/2 + 1;
        endIdx = nAreas;                        
    end            
    lesionCntr = centres(lesionIdx,:);
    [D,I] = pdist2(centres(startIdx:endIdx,:),lesionCntr,'euclidean','smallest',nNeighbors);
    I = (startIdx-1) + I;
    lSC(I,:) = 0;
    lSC(:,I) = 0;
    
%     [pathstr,subj_name,ext] = fileparts(scPath);
%     t = strsplit(subj_name,'_');
%     subj_foldr = t{1};
%     clear t;
    savePath = [saveDir,'sc_l',num2str(lesionIdx),'_nn',num2str(nNeighbors),'.hdf5'];
%     savePath = [saveDir,'sc_',num2str(lesionIdx),'_',strtrim(deblank(regionNames{lesionIdx})),'_nn5.hdf5'];
%     saveFigPath = [saveDir,'sc_',num2str(lesionIdx),'_',strtrim(deblank(regionNames{lesionIdx})),'_nn5.jpeg'];
%     savePath = [saveDir,subj_foldr,'/','sc_',num2str(lesionIdx),'_',strtrim(deblank(regionNames{lesionIdx})),'_nn',num2str(nNeighbors),'.hdf5'];    
    rmCmd = ['rm ',savePath];
    if(~system(['test -e ' savePath]))
        system(rmCmd);
    end
    fprintf('Saving Lesioned SC: %s\n',savePath);
    h5create(savePath,'/C',[nAreas nAreas]);
    h5write(savePath,'/C',lSC);
    h5create(savePath,'/Neighbours',length(I));
    h5write(savePath,'/Neighbours',I);
    h5writeatt(savePath,'/C','LesionName',lesionName);
    h5writeatt(savePath,'/C','LesionIdx',lesionIdx);
%     fprintf('Plotting Lesioned SC: %s\n\n',saveFigPath);
%     figure('Visible','off');
%     plot_lesioned_sc(savePath,'LesionColorMap_998.mat',regionNames);
%     print('-djpeg',saveFigPath);
%     close all;
end