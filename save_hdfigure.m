function [] = save_hdfigure(savePath)

%%%%%%%% plot the figure %%%%%%%
% imagesc(sc_lpcun)
% colormap(clrMap)
% colorbar()
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%


width = 3;     % 3  % Width in inches
height = 2;  % 2  % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
pos = get(gcf, 'Position');
set(gcf,'Units','inches')
set(gcf, 'Position', [pos(1) pos(2) width height]); %<- Set size
set(gca, 'LineWidth', alw); %<- Set properties
% set(gca,'Units','normalized','Position',[0.18 0.2 0.6 0.6])
% cpos = get(c,'Position');
% set(c,'Units','normalized','Position',[0.85 0.2 cpos(3) 0.6])

% Set Tick Marks
% set(gca,'XTick',[24],'YTick',[24]);
% xtklbls = cell(1);
% xtklbls{1} = 'PCUN';
% ytklbls = cell(1);
% ytklbls{1} = 'PCUN';
% set(gca,'XTickLabel',ytklbls,'YTickLabel',ytklbls)

set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
% set(gcf,'PaperPositionMode', 'auto');


% Save the file as PNG
print([savePath '.png'],'-dpng','-r300');
savefig([savePath '.fig']);
disp(['saved: ' savePath '.png'])
end