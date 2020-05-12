DMNRegions = [13 22 7 28 27 8 15 47 41 62 61 42 49];
nDMNRegions = length(DMNRegions);
scPaths = {'../../data/connectivity_68/lesioned/sc_24_rPCUN.hdf5',...
           '../../data/connectivity_68/lesioned/sc_2_rCAC.hdf5',...
           '../../data/connectivity_68/lesioned/sc_7_rIP.hdf5',...
           '../../data/connectivity_68/lesioned/sc_30_rSMAR.hdf5',...
           '../../data/connectivity_68/lesioned/sc_29_rST.hdf5',...
           '../../data/connectivity_68/lesioned/sc_8_rIT.hdf5',...
           '../../data/connectivity_68/lesioned/sc_27_rSF.hdf5',...
           '../../data/connectivity_68/lesioned/sc_3_rCMF.hdf5',...
           '../../data/connectivity_68/lesioned/sc_17_rPOPE.hdf5',...
           '../../data/connectivity_68/lesioned/sc_10_rLOCC.hdf5',...
           '../../data/connectivity_68/lesioned/sc_23_rPREC.hdf5',...
           '../../data/connectivity_68/lesioned/sc_58_lPCUN.hdf5',...
           '../../data/connectivity_68/lesioned/sc_36_lCAC.hdf5',...
           '../../data/connectivity_68/lesioned/sc_41_lIP.hdf5',...
           '../../data/connectivity_68/lesioned/sc_64_lSMAR.hdf5',...
           '../../data/connectivity_68/lesioned/sc_63_lST.hdf5',...
           '../../data/connectivity_68/lesioned/sc_42_lIT.hdf5',...
           '../../data/connectivity_68/lesioned/sc_61_lSF.hdf5',...
           '../../data/connectivity_68/lesioned/sc_37_lCMF.hdf5',...
           '../../data/connectivity_68/lesioned/sc_51_lPOPE.hdf5',...
           '../../data/connectivity_68/lesioned/sc_44_lLOCC.hdf5',...
           '../../data/connectivity_68/lesioned/sc_57_lPREC.hdf5'};
for i = 1:length(scPaths)
    nbrs = h5read(scPaths{i},'/Neighbours');
    DMNExtent = length(intersect(nbrs,DMNRegions))/nDMNRegions;
    disp(strcat(['DMN extent for ',scPaths{i},' :',num2str(DMNExtent)]));
end