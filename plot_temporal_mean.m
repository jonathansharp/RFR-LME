% Plot temporal mean of spatial surface variable
function plot_temporal_mean(lim,dim,lat,lon,z,zmin,zmax,lev,varname,lab,reg)

figure('visible','off'); worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
pcolorm(repmat(lat',dim.x,1),repmat(lon,1,dim.y),mean(z,3,'omitnan'));
land = shaperead('landareas', 'UseGeoCoords',true);
% geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(parula(lev));
caxis([zmin,zmax]);
c.TickLength = 0;
c.Label.String = lab;
cbarrow;

% save figure
exportgraphics(gcf,['Figures/' reg '_' varname '.png']);
close