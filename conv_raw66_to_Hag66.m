function [ scHag ] = conv_raw66_to_Hag66( scRaw, raw66Hag66mapPath )
% converts a raw order SC into Hagmann order SC
    nAreas = size(scRaw,1);
    scHag = zeros(nAreas,nAreas);
    raw66Hag66map = dlmread(raw66Hag66mapPath);
    for i = 1:nAreas
        for j = 1:nAreas
            scHag(raw66Hag66map(i)+1,raw66Hag66map(j)+1) = scRaw(i,j);
        end
    end
end

