load neuro_act.mat;
clearvars -except wEI;
scPath = '../../data/connectivity_68/lesioned/sc_24_rPCUN_nn1.hdf5';
simTime = 8*60*1000;
dt = 1;
G = 0.6;
noise_amp = 0.001;
lesioned_areas = [24];
locc_stims = 0.02:0.01:0.1;
nstims = length(locc_stims);
nAreas = 68;
avgpwr_prior = zeros(nstims,nAreas);
alpha_pwr_prior = zeros(nstims,nAreas);
beta_pwr_prior = zeros(nstims,nAreas);
gamma_pwr_prior = zeros(nstims,nAreas);
avgpwr_post = zeros(nstims,nAreas);
alpha_pwr_post = zeros(nstims,nAreas);
beta_pwr_post = zeros(nstims,nAreas);
gamma_pwr_post = zeros(nstims,nAreas);
for i = 1:length(locc_stims)
    Iext = zeros(68,1);
    Iext([10 44]) = locc_stims(i);
    [fc_prior, maxfrN_prior, curr2_prior, neuro_act_prior, wEI_prior]=DMF_dev(scPath,'/C',false,false,simTime,dt,G,noise_amp,lesioned_areas,wEI,true,Iext);
    [fc_post, maxfrN_post, curr2_post, neuro_act_post, wEI_post]=DMF_dev(scPath,'/C',true,false,simTime,dt,G,noise_amp,lesioned_areas,wEI,true,Iext);
    % [fc_post2, maxfrN_post2, curr2_post2, neuro_act_post2, wEI_post2]=DMF_dev(scPath,'/C',false,false,simTime,dt,G,noise_amp,lesioned_areas,wEI_post1,true,Iext);
    for j = 1:nAreas
        avgpwr_prior(i,j) = bandpower(neuro_act_prior(5000:end,j));
        % avgpwr_post1 = bandpower(neuro_act_post1(5000:end,1));
        avgpwr_post(i,j) = bandpower(neuro_act_post(5000:end,j));
        alpha_pwr_prior(i,j) = bandpower(neuro_act_prior(5000:end,j),1000,[7 15]);
        alpha_pwr_post(i,j) = bandpower(neuro_act_post(5000:end,j),1000,[7 15]);
        beta_pwr_prior(i,j) = bandpower(neuro_act_prior(5000:end,j),1000,[16 31]);
        beta_pwr_post(i,j) = bandpower(neuro_act_post(5000:end,j),1000,[16 31]);
        gamma_pwr_prior(i,j) = bandpower(neuro_act_prior(5000:end,j),1000,[32 100]);
        gamma_pwr_post(i,j) = bandpower(neuro_act_post(5000:end,j),1000,[32 100]);    
    end
end
save('../../results/68/power_freqbands.mat');