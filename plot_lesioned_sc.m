function [c] = plot_lesioned_sc(scPath,clrMapPath,regionNames)
    sc = h5read(scPath,'/C');
    nbrs = h5read(scPath,'/Neighbours');
%     nAreas = size(sc,1);
    clrMap = load(clrMapPath);
    clrMap = clrMap.clrMap;
%     lesionIdx = h5readatt(scPath,'/C','LesionIdx');
%     sc(nbrs,:) = -0.5;
%     sc(:,nbrs) = -0.5;
    sc(nbrs,:) = -0.05;
    sc(:,nbrs) = -0.05;
    imagesc(sc,[-0.05 max(max(sc))]);
    colormap(clrMap);
    maxVal = max(max(sc));
    c = colorbar('YTick',[-0.025,0,maxVal/4,maxVal/2,(maxVal/2+maxVal)/2,maxVal],'YTickLabel',{'Lesion','0',num2str(maxVal/4,'%.2f'),num2str(maxVal/2,'%.2f'),num2str((maxVal/2+maxVal)/2,'%.2f'),num2str(maxVal,'%.2f')});
    
%     colorbar('YTick',[-0.15,0,0.3,0.6,0.8,1],'YTickLabel',{'Lesion','0','0.2','0.4','0.6','0.8','1'});
%     nbrs = sort(nbrs);
%     xTickLabel = cell(length(nbrs),1);
%     yTickLabel = cell(length(nbrs),1);
%     for i = 1:length(nbrs)
% %         if(lesionIdx <= nAreas/2)
% %             xTickLabel{i} = strcat(['r',regionNames{nbrs(i)}]);
% %             yTickLabel{i} = strcat(['r',regionNames{nbrs(i)}]);
% %         else
% %             xTickLabel{i} = strcat(['l',regionNames{nbrs(i)}]);
% %             yTickLabel{i} = strcat(['l',regionNames{nbrs(i)}]);
% %         end
%         xTickLabel{i} = strcat(regionNames{nbrs(i)});
%         yTickLabel{i} = strcat(regionNames{nbrs(i)});
%     end
%     set(gca,'XTick',nbrs,'XTickLabel',xTickLabel);
%     set(gca,'YTick',nbrs,'YTickLabel',yTickLabel);
%     set(gca,'XTick',lesionIdx,'XTickLabel',regionNames{lesionIdx});
%     set(gca,'YTick',lesionIdx,'YTickLabel',regionNames{lesionIdx});
    axpos = get(gca,'Position');
    cpos = get(c,'Position');
    cpos(3) = 0.5*cpos(3);
    set(c,'Position',cpos);
    set(gca,'Position',axpos);
end