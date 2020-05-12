function [] = gen_diffuse_lesioned_sc(scPath,saveDir)
    sc = h5read(scPath,'/C');
    nAreas = size(sc,1);
    nConns = nAreas*nAreas;
    for i = 1:50
        t = sc;
        randIdxs = randperm(nConns,floor((i/100)*nConns));
        t(randIdxs) = 0;
        fName = [saveDir,'/','sc_diffuse_',int2str(i),'.hdf5'];
        if(exist(fName,'file'))
            system(['rm ',fName]);
        end
        h5create(fName,'/C',size(sc));
        h5write(fName,'/C',t);
        h5create(fName,'/Lesioned_Idxs',size(randIdxs));
        h5write(fName,'/Lesioned_Idxs',randIdxs);
    end
end