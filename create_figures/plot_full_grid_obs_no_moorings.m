%% Plot the percentage of grid cells with data
load('Data/socat_gridded_2023_no_moorings','SOCAT_grid');
figure('visible','off');
worldmap([SOCAT_grid.lim.latmin ...
    SOCAT_grid.lim.latmax],...
   [SOCAT_grid.lim.lonmin ...
    SOCAT_grid.lim.lonmax]);
setm(gca,'MapProjection','miller');
set(gca,'fontsize',16);
pcolorm(repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1),...
        repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y),...
        SOCAT_grid.num_months);
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
mycolormap = jet(16);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-1 31]);
c.TickLength = 0;
c.Label.String = 'Number of Months Represented';
cbarrow('up');
mlabel off;
plabel off;
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/obs_no_moorings.png');
close
% clean up
clear