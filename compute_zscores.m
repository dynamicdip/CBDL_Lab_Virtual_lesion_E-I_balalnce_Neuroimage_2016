function fc_zscores = compute_zscores(matA,matB)
    df = 230;
    nAreas = size(matA,1);
    matA(1:nAreas+1:nAreas*nAreas) = 0;
    matB(1:nAreas+1:nAreas*nAreas) = 0;
    matA = atanh(matA);
    matB = atanh(matB);
    fc_zscores = (matA - matB) ./ (sqrt(2/(df - 3)));
end