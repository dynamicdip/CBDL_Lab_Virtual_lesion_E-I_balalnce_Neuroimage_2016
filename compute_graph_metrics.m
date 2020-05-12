function compute_graph_metrics(connectivityType,connectivityPath,hdf5SavePath)
    addpath /home/iiit/anirudh/BCT/2015_01_25' BCT';
%     regionNamesPath = '../../data/connectivity_68/fs_region_abbrev_sort.txt';
%     scHCPath = '../../data/connectivity_68/Human_68.hdf5';
%     lscPath = '../../data/connectivity_68/lesioned/sc_24_rPCUN.hdf5';
%     fcHCPath = '../../data/connectivity_68/Human_68.hdf5';
%     lfcwoficPath = '../../results/68/fc_24_rPCUN_G0.15_wofic_nn4.hdf5';
%     lfcficPath = '../../results/68/fc_24_rPCUN_G0.5_fic_nn4.hdf5';
%     regionNames = strsplit(fileread(regionNamesPath),'\n');
    
%     lsc = h5read(lscPath,'/C');
%     scHClen = weight_conversion(scHC, 'lengths');
%     lsclen = weight_conversion(lsc, 'lengths');
    rmCmd = strcat(['rm ',hdf5SavePath]);
    if(~system(strcat(['test -e ',hdf5SavePath])))
        system(rmCmd);
    end
    if(strcmp(connectivityType,'structural'))
        sc = h5read(connectivityPath,'/C');
        scLen = weight_conversion(sc, 'lengths');
        %Clustering Coefficient
        clstrCoef = clustering_coef_wd(sc);
        h5create(hdf5SavePath,'/ClusteringCoefficient',size(clstrCoef));
        h5write(hdf5SavePath,'/ClusteringCoefficient',clstrCoef);
        % subplot(1,2,1)
        % [clstrCoefHCsort clstrCoefHCidxs] = sort(clstrCoefHC);
        % figure(1);
        % subplot(2,1,1);
        % bar(clstrCoefHCsort);
%         clstrCoeflsc = clustering_coef_wd(lsc);

        % Characteristic Path length and Global Efficiency
        D = distance_wei(scLen);
        [lambda,efficiency] = charpath(D);
        h5create(hdf5SavePath,'/PathLength',size(lambda));
        h5write(hdf5SavePath,'/PathLength',lambda);
        h5create(hdf5SavePath,'/GlobalEfficiency',size(efficiency));
        h5write(hdf5SavePath,'/GlobalEfficiency',efficiency);
        % Dlsc = distance_wei(lsc);
%         Dlsc = distance_wei(lsclen);
%         [lambdalsc,efficiencylsc] = charpath(Dlsc);

        % Degree
        [ideg,odeg,deg] = degrees_dir(sc);
        h5create(hdf5SavePath,'/InDegree',size(ideg));
        h5write(hdf5SavePath,'/InDegree',ideg);
        h5create(hdf5SavePath,'/OutDegree',size(odeg));
        h5write(hdf5SavePath,'/OutDegree',odeg);
        h5create(hdf5SavePath,'/Degree',size(deg));
        h5write(hdf5SavePath,'/Degree',deg);
%         [ideglsc,odeglsc,deglsc] = degrees_dir(lsc);

        % Assortativity
        r = assortativity_wei(sc,3);
        h5create(hdf5SavePath,'/Assortativity',size(r));
        h5write(hdf5SavePath,'/Assortativity',r);
%         rlsc = assortativity_wei(lsc,3);

        % Node Betweenness Centrality
        nbc = betweenness_wei(scLen);
        h5create(hdf5SavePath,'/NodeBetweennessCentrality',size(nbc));
        h5write(hdf5SavePath,'/NodeBetweennessCentrality',nbc);
%         nbclsc = betweenness_wei(lsclen);

        % edge Betweeness Centrality
        ebc = edge_betweenness_wei(scLen);
        h5create(hdf5SavePath,'/EdgeBetweennessCentrality',size(ebc));
        h5write(hdf5SavePath,'/EdgeBetweennessCentrality',ebc);
%         ebclsc = edge_betweenness_wei(lsclen);

        % Modularity
        [m,q] = community_louvain(sc);
        h5create(hdf5SavePath,'/CommunityStructure',size(m));
        h5write(hdf5SavePath,'/CommunityStructure',m);
        h5create(hdf5SavePath,'/Modularity',size(q));
        h5write(hdf5SavePath,'/Modularity',q);
%         [mlsc,qlsc] = community_louvain(lsc);

        % Connected Components
        scBin = weight_conversion(sc, 'binarize');
        [scComps,scCompSizes] = get_components(scBin);
        h5create(hdf5SavePath,'/ConnectedComponents',size(scComps));
        h5write(hdf5SavePath,'/ConnectedComponents',scComps);
        h5create(hdf5SavePath,'/ConnectedComponentSizes',size(scCompSizes));
        h5write(hdf5SavePath,'/ConnectedComponentSizes',scCompSizes);
%         lscBin = weight_conversion(lsc, 'binarize');
%         [lscComps,lscCompSizes] = get_components(lscBin);

        % Comparison of network community partitions
%         [VInscHCvssclsc, MInscHCvssclsc]       = partition_distance(mHC,mlsc);

        % Rich Club Coefficient
        [Rw] = rich_club_wd(sc);
        h5create(hdf5SavePath,'/RichClubCoefficient',size(Rw));
        h5write(hdf5SavePath,'/RichClubCoefficient',Rw);
    elseif(strcmp(connectivityType,'functional'))
        %%%%%%%%%Functional Connectivity%%%%%%%%%%%%%%%%
        fc = h5read(connectivityPath,'/CC');
%         fcHC = h5read(fcHCPath,'/CC');
%         lfcwofic = h5read(lfcwoficPath,'/avgCC');
%         lfcfic = h5read(lfcficPath,'/avgCC');

        % Modularity
        [Ci, Q] = modularity_louvain_und_sign(fc,'sta');
        h5create(hdf5SavePath,'/CommunityStructure',size(Ci));
        h5write(hdf5SavePath,'/CommunityStructure',Ci);
        h5create(hdf5SavePath,'/Modularity',size(Q));
        h5write(hdf5SavePath,'/Modularity',Q);
%         [Cilfcwofic, Qlfcwofic] = modularity_louvain_und_sign(lfcwofic,'sta');
%         [Cilfcfic, Qlfcfic] = modularity_louvain_und_sign(lfcfic,'sta');

        % Nodal Diversity Coefficient
%         [Hpos, Hneg] = diversity_coef_sign(fc,Ci);        
%         [Hposlfcwofic, Hneglfcwofic] = diversity_coef_sign(lfcwofic,Cilfcwofic);
%         [Hposlfcfic, Hneglfcfic] = diversity_coef_sign(lfcfic,Cilfcfic);


        % Nodal Participation Coefficient
%         [Ppos, Pneg] = participation_coef_sign(fc,Ci);
%         [Pposlfcwofic, Pneglfcwofic] = participation_coef_sign(lfcwofic,Cilfcwofic);
%         [Pposlfcfic, Pneglfcfic] = participation_coef_sign(lfcfic,Cilfcfic);

        % Partition distance
%         [VInHCvslfcwofic MInHCvslfcwofic] = partition_distance(CiHC,Cilfcwofic);
%         [VInHCvslfcfic MInHCvslfcfic] = partition_distance(CiHC,Cilfcfic);

        % connected components
        fcBin = weight_conversion(fc, 'binarize');
        [fcComps,fcCompSizes] = get_components(fcBin);
        h5create(hdf5SavePath,'/ConnectedComponents',size(fcComps));
        h5write(hdf5SavePath,'/ConnectedComponents',fcComps);
        h5create(hdf5SavePath,'/ConnectedComponentSizes',size(fcCompSizes));
        h5write(hdf5SavePath,'/ConnectedComponentSizes',fcCompSizes);
%         lfcBin = weight_conversion(lfcwofic, 'binarize');
%         [lfcComps,lfcCompSizes] = get_components(lfcBin);
        
        % Rich Club Coefficient
        [Rw] = rich_club_wu(fc);
        h5create(hdf5SavePath,'/RichClubCoefficient',size(Rw));
        h5write(hdf5SavePath,'/RichClubCoefficient',Rw);
    end
    
end

