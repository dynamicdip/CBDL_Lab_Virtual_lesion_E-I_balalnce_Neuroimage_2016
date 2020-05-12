function map_998_to_66(scPath998,scName,savePath66)
    sc_998 = h5read(scPath998,scName);
    lbls_map = h5read('../../data/connectivity_998/Hagmann_998.hdf5','/Map_998_to_66');
    sc_66 = zeros(66,66);
    for i = 1:66
        trgt_regn_idxs = find(lbls_map == i);
        for j = 1:66
            src_regn_idxs = find(lbls_map == j);
%             sc_66(i,j) = sum(sum(sc_998(src_regn_idxs,trgt_regn_idxs))) / length(trgt_regn_idxs);
            sc_66(i,j) =  mean(mean(sc_998(src_regn_idxs,trgt_regn_idxs)));
        end
    end
    hdf5write(savePath66,'/C_resampled_roi_avg',sc_66);
end


% hag_data = load('../../data/connectivity_998/DSI_enhanced.mat');
% sc_998_fbden = hag_data.CIJ_fbden_average;
% sc_998_resampled = hag_data.CIJ_resampled_average;
% fc_998 = hag_data.COR_fMRI_average;
% sc_66_fbden = zeros(66,66);
% sc_66_resampled = zeros(66,66);
% fc_66 = zeros(66,66);
% roi_talairach = zeros(66,3);
% for i = 1:66
%     trgt_regn_idxs = find(hag_data.roi_lbls == i);
%     for j = 1:66
%         src_regn_idxs = find(hag_data.roi_lbls == j);
%         sc_66_fbden(i,j) = mean(mean(sc_998_fbden(src_regn_idxs,trgt_regn_idxs)));
%         sc_66_resampled(i,j) =  mean(mean(sc_998_resampled(src_regn_idxs,trgt_regn_idxs)));
% %         if(i<j)
% %             fc_66(i,j) = mean(mean(fc_998(src_regn_idxs,trgt_regn_idxs)));
% %         elseif(i == j)
% %             fc_66(i,j) = 1;
% %         else
% %             fc_66(i,j) = fc_66(j,i);
% %         end
%         fc_66(i,j) = mean(mean(fc_998(src_regn_idxs,trgt_regn_idxs)));
%     end
%     roi_talairach(i,:) = mean(hag_data.talairach(trgt_regn_idxs,:));
%     roi_lbls{i} = hag_data.anat_lbls(i,:);
% end
% % save('../../data/connectivity_66/Hagmann_66.mat','sc_66_fbden','sc_66','fc_66','roi_lbls','roi_talairach');
% hdf5write('../../data/connectivity_66/Hagmann_66.hdf5','/C_fbden_roi_avg',sc_66_fbden);
% hdf5write('../../data/connectivity_66/Hagmann_66.hdf5','/C_resample_roi_avg',sc_66_resampled,'WriteMode','append');
% hdf5write('../../data/connectivity_66/Hagmann_66.hdf5','/CC',fc_66,'WriteMode','append');
% hdf5write('../../data/connectivity_66/Hagmann_66.hdf5','/Talairach',roi_talairach,'WriteMode','append');
% hdf5write('../../data/connectivity_66/Hagmann_66.hdf5','/ROI_lbls',roi_lbls,'WriteMode','append')
