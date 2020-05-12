function save_hdf5(savePath,dsName,data)
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create(savePath,'H5F_ACC_TRUNC',fcpl,fapl);
    type_id = H5T.copy('H5T_NATIVE_DOUBLE');
    dims = size(data);
    h5_dims = fliplr(dims);
    h5_maxdims = h5_dims;
    data_rank = length(size(data));
    space_id = H5S.create_simple(data_rank,h5_dims,h5_maxdims);
    lcpl = H5P.create('H5P_LINK_CREATE');
    dcpl = H5P.create('H5P_DATASET_CREATE');
    dapl = H5P.create('H5P_DATASET_ACCESS');
    dset_id = H5D.create(fid,dsName,type_id,space_id,lcpl,dcpl,dapl);
    file_space_id = H5D.get_space(dset_id);
    H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data);
    %     H5D.write(dset_id,type_id,space_id,fid,data);
    H5S.close(space_id);
    H5T.close(type_id);
%     H5D.close(file_space_id);
    H5D.close(dset_id);
    H5F.close(fid);
    H5P.close(dcpl);
    H5P.close(fcpl);
    H5P.close(fapl);
    
end