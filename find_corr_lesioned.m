function [ xyCorr ] = find_corr_lesioned(x, y, nbrs)
    % computes the correlation between two matrices by finding the
    % correlation of lower diagonal elements of each matrix
    nAreas = size(x,1);
    if(nAreas ~= size(y,1))
       disp('Both the matrices must be of same size') 
       return
    end
    IsubDiag = find(tril(ones(nAreas),-1)); % indices of lower diagonal elements of SC matrix
    tmp = zeros(nAreas,nAreas);
    tmp(nbrs,:) = 1;
    tmp(:,nbrs) = 1;
    nbrsIdxs = find(tmp == 1);
%     for i = 1:length(nbrs)
%         IsubDiag = IsubDiag((mod((IsubDiag - nbrs(i)),nAreas) ~= 0));
%         IsubDiag = IsubDiag(~((IsubDiag-(nbrs(i)-1)*nAreas)>=1 & (IsubDiag-(nbrs(i)-1)*nAreas)<=nAreas));
%     end
    IsubDiag = setdiff(IsubDiag,nbrsIdxs);
    xyCorr=corrcoef(atanh(x(IsubDiag)),atanh(y(IsubDiag)));
    xyCorr = xyCorr(1,2);
end