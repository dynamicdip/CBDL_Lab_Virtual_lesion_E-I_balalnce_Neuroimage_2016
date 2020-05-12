scPath = '../../data/connectivity_68/Human_68.hdf5'
sc = h5read(scPath,'/C');
fic = true
ffi = false
simTime = 8*60*1000
dt = 0.1
G = 0.6
noiseAmp = 0.001
lesionAreas = []
ficWeights = []
setFICWeights = false
p = parpool(12);
fcs = cell(10);
frNs = cell(10);
parfor i = 1:10
    [fcs{i},frNs{i}]=DMF_lesioned(scPath,fic,ffi,simTime,dt,G,noiseAmp,lesionAreas,ficWeights,setFICWeights);
    frNs{i} = max(frNs{i});
end
%fc = zeros(size(sc));
% for i = 1:5
%     fc = fc + fcs{i};
% end
%fc = fc ./ 5;

save('../../results/68/sim_healthy_1.mat','fcs','frNs')
p.delete;
clear p;
