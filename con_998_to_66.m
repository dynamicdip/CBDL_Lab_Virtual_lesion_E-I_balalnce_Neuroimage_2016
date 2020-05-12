function [ SC_66 ] = con_998_to_66( scPath,regionIdxsPath, raw66Hag66mapPath )
% converts the Strctural Connectivity of 998 areas to 66 areas as mentioned
% in: How local excitation-inhibition ratio impacts whole brain dynamics Deco et.al
% Note: the order of areas are according to Hagmann order
    SC_998 = h5read(scPath,'/C');
    SC_66 = zeros(66,66);
    regionIdxs = dlmread(regionIdxsPath);
    for i = 1:66
            if(i ~= 66)
                nSourceROIs = (regionIdxs(i+1) - regionIdxs(i));
            else
                nSourceROIs = 998 - regionIdxs(i);
            end
            rowStart = regionIdxs(i) + 1;
            rowEnd = regionIdxs(i) + nSourceROIs;
        for j = 1:66
            if(i == j)
                continue; 
            else
                if(j ~= 66)
                    nTargetROIs = regionIdxs(j+1) - regionIdxs(j);
                else
                    nTargetROIs = 998 - regionIdxs(j);
                end
                colStart = regionIdxs(j) + 1;
                colEnd = regionIdxs(j) + nTargetROIs;
                SC_66(i,j) = sum(sum(SC_998(rowStart:rowEnd,colStart:colEnd))) / nTargetROIs;
%                 SC_66(i,j) = sum(sum(SC_998(rowStart:rowEnd,colStart:colEnd)));
            end
        end
    end
    SC_66 = SC_66';
    SC_66 = conv_raw66_to_Hag66(SC_66, raw66Hag66mapPath);
end

