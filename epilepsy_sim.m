clear all;
sc = h5read('../../data/connectivity_68/Human_66.hdf5','/C');
roiNos = [13,47,24,58,22,56];
p = parpool(12);
sf = 0:0.1:1;
fcs = cell(length(sf));
frN = cell(length(sf));
scScaled = cell(length(sf));
parfor i = 1:length(sf)
    scScaled{i} = scale_sc_rois(sc,roiNos,sf(i));
    [fcs{i}, frN{i}] = DMF_epilepsy(scScaled{i},true,false,480000,1,0.6,0.001);
end
save('fcs_epilepsy.mat','fcs');