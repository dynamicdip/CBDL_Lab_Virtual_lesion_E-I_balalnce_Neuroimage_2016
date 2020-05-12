function [ xyCorr ] = find_corr(x, y)
    % computes the correlation between two matrices by finding the
    % correlation of lower diagonal elements of each matrix
    nAreas = size(x,1);
    if(nAreas ~= size(y,1))
       disp('Both the matrices must be of same size') 
       return
    end
    IsubDiag = find(tril(ones(nAreas),-1)); % indices of lower diagonal elements of SC matrix
    xyCorr=corrcoef(atanh(x(IsubDiag)),atanh(y(IsubDiag)));
    xyCorr = xyCorr(1,2);
end