function lesionStrength = lesion_strength(sc,lesionAreas)
%% computes the strength of the lesion which is the sum of weights of all connections
%% to lesioned nodes
lesionStrength = sum(sum(sc(lesionAreas,:)));
%lesionStrength = lesionStrength + sum(sum(sc(:,lesionAreas)));

end