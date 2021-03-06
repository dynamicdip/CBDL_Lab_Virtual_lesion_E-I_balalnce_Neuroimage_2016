clear all;
scPath = '../../data/connectivity_114/cmp_114.hdf5';
hdf5savePath = '../../data/connectivity_114/cmp_114.hdf5';
matSavePath = '../../data/connectivity_114/cmp_114.mat';
mu = 0.5;
sigma = 0.1;
sc = h5read(scPath,'/FiberCount');
nonzero_idxs = find(sc ~= 0);
scvals = sc(nonzero_idxs);
nAreas = size(sc,1);
normal_samples = sigma*randn(1,length(scvals)) + mu;
sc_resampled = zeros(nAreas,nAreas);
[t1,idx1] = sort(scvals);
[t2,idx2] = sort(normal_samples);
sc_resampled(nonzero_idxs(idx1)) = normal_samples(idx2);
hdf5write(hdf5savePath,'/C_resample',sc_resampled,'WriteMode','append');
save(matSavePath,'sc_resampled','-append');