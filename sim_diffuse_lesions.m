scPath = '../../data/connectivity_68/lesioned/diffuse/sc_diffuse_'
fcPath = '../../results/68/diffuse/fc_post_diffuse_'
p = parpool(6);
fcs = cell(5);
frNs = cell(5);
nAreas = 68
nItersPerSim = 5
simTime = 8*60*1000
fic = true
ffi = false
dt = 0.1
G = 0.6
noiseAmp = 0.001
t = load('sim_hc.mat','wEI');
ficWeights = t.wEI;
clear t;
fc_emp = h5read('../../data/connectivity_68/Human_68.hdf5','/CC');
for i = 1:50
    fcCorrs = zeros(nItersPerSim,1);
    parfor j = 1:nItersPerSim
        [fcs{j},frNs{j}] = DMF_lesioned([scPath,int2str(i),'.hdf5'],fic,ffi,simTime,dt,G,noiseAmp,[],ficWeights,true);
        fcCorrs(j) = find_corr(fc_emp,fcs{j});
    end
    maxfrNs = zeros(nItersPerSim,nAreas);
    maxFCIdx = 0;
    for k = 1:nItersPerSim
        if (k == 1)
           avgFC = fcs{k};
           maxFCIdx = k;
        else
           avgFC = avgFC + fcs{k};
           if(fcCorrs(k) < fcCorrs(maxFCIdx))
              maxFCIdx = k; 
           end
        end
        maxfrNs(k,:) = max(frNs{k});
    end
    avgFC = avgFC ./ nItersPerSim;
    if(exist([fcPath,int2str(i),'.hdf5'],'file'))
       system(['rm ',[fcPath,int2str(i),'.hdf5']]);
    end
    h5create([fcPath,int2str(i),'.hdf5'],'/avgCC',size(avgFC));
    h5write([fcPath,int2str(i),'.hdf5'],'/avgCC',avgFC);
    h5create([fcPath,int2str(i),'.hdf5'],'/maxCC',size(fcs{maxFCIdx}));
    h5write([fcPath,int2str(i),'.hdf5'],'/maxCC',fcs{maxFCIdx});
    h5writeatt([fcPath,int2str(i),'.hdf5'],'/maxCC','maxFCIdx',maxFCIdx);
    h5create([fcPath,int2str(i),'.hdf5'],'/maxfrNs',size(maxfrNs));
    h5write([fcPath,int2str(i),'.hdf5'],'/maxfrNs',maxfrNs);
end