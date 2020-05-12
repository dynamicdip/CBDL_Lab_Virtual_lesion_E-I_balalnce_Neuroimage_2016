% Create Normalized Structural Connectivity matrix for 68 areas in Desikan
% Killiany Atlas
load('../../data/connectivity_114/sc_cmp_scale60.mat');
load('../../data/connectivity_114/surfareas.mat');
start_idx_rh = 1;
end_idx_rh = 57;
start_idx_lh = 65;
end_idx_lh = 121;
for i = start_idx_rh:end_idx_rh
    cmp.node_names{i} = nodes.dn_name{i};
end
for i = start_idx_lh:end_idx_lh
    cmp.node_names{i-(start_idx_lh-end_idx_rh-1)} = nodes.dn_name{i};
end
cmp.number_of_fibers = sc.number_of_fibers([start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh],[start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh]);
cmp.fiber_length_mean = sc.fiber_length_mean([start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh],[start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh]);
% zero_idxs=find(cmp_68.fiber_length_mean == 0);
nAreas = size(cmp.number_of_fibers,1);
cmp.sc_norm = zeros(nAreas,nAreas);
cmp.surf_areas = zeros(1,nAreas);
for i = 1:nAreas
    for j = 1:nAreas
        if(strcmp(cmp.node_names{i},surfareas{j,1}))
            cmp.surf_areas(i) = surfareas{j,2};
            break;
        end
    end
end
for i = 1:nAreas
    for j = 1:nAreas
        if(cmp.fiber_length_mean(i,j) ~= 0)
            cmp.sc_norm(i,j) = (2 / (cmp.surf_areas(i) + cmp.surf_areas(j))) * cmp.number_of_fibers(i,j) * (1/cmp.fiber_length_mean(i,j));
        else
            cmp.sc_norm(i,j) = (2 / (cmp.surf_areas(i) + cmp.surf_areas(j))) * cmp.number_of_fibers(i,j);
        end
        cmp.sc_norm1(i,j) = (2 / (cmp.surf_areas(i) + cmp.surf_areas(j))) * cmp.number_of_fibers(i,j);
    end
end
cmp.sc_norm1 = cmp.sc_norm1 ./ max(max(cmp.sc_norm1));

% Create FC matrix for 68 areas in Desikan Killinay atlas
fc_cmp = load('../../data/connectivity_114/fc_cmp_scale60.mat');
cmp.fc = fc_cmp.sc.corr([start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh],[start_idx_rh:end_idx_rh start_idx_lh:end_idx_lh]);

save('../../data/connectivity_114/cmp_114.mat','cmp');
hdf5write('../../data/connectivity_114/cmp_114.hdf5','/C',cmp.sc_norm,'/C1',cmp.sc_norm1,'/FiberCount',cmp.number_of_fibers,...
          '/MeanFiberLength',cmp.fiber_length_mean,'/SurfAreas',cmp.surf_areas,'/CC',cmp.fc);
