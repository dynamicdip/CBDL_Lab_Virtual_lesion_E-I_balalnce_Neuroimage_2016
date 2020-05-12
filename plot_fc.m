function plot_fc(fc)
%     fc = h5read(fcPath,'/cc');
    imagesc(fc);
    colormap(darkb2r(-0.5,1));
    colorbar();
    set(gca,'Fontsize', 20);
end