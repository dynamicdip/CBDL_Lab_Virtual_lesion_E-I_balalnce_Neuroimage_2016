fcPath = '../../results/68/diffuse/fc'
simHC = load('sim_hc.mat');
fcHCSim = simHC.fc;
clear simHC;
fcEmp = h5read('../../data/connectivity_68/Human_68.hdf5','/CC');
fcDistsPrior = zeros(1,50);
fcDistsPost = zeros(1,50);
for i = 1:50
    fcPrior = h5read([fcPath,'_prior_diffuse_',int2str(i),'.hdf5'],'/avgCC');
    fcPost = h5read([fcPath,'_post_diffuse_',int2str(i),'.hdf5'],'/avgCC');
    fcDistsPrior(i) = fc_distance(fcEmp,fcPrior,[]);
    fcDistsPost(i) = fc_distance(fcEmp,fcPost,[]);
end
plot(1:50,fcDistsPrior,'ko','MarkerFaceColor','k');
hold on;
plot(1:50,fcDistsPost,'bo','MarkerFaceColor','b');