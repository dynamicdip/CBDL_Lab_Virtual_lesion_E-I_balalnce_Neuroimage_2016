function [scNew] = scale_sc_rois(scOld,roiNos,scale)

scNew = scOld;
scNew(:,roiNos) = scale*scNew(:,roiNos);
scNew(roiNos,:) = scale*scNew(roiNos,:);

end