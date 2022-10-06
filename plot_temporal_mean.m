% Plot temporal mean of spatial surface variable
% 
% Written by J.D. Sharp: 8/24/22
% Last updated by J.D. Sharp: 9/15/22
% 

function plot_temporal_mean(lim,dim,lat,lon,z,zmin,zmax,lev,varname,lab,reg)

figure('visible','off'); worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
pcolorm(repmat(lat',dim.x,1),repmat(lon,1,dim.y),mean(z,3,'omitnan'));
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(parula(lev));
caxis([zmin,zmax]);
c.TickLength = 0;
c.Label.String = lab;
cbarrow;

% save figure
if ~isfolder(['Figures/' reg]); mkdir(['Figures/' reg]); end
exportgraphics(gcf,['Figures/' reg '/' varname '.png']);
close