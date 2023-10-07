%% Plot detrended gridded mean pCO2
load('Data/socat_gridded_2023','SOCAT_grid');
figure('visible','off');
worldmap([SOCAT_grid.lim.latmin ...
    SOCAT_grid.lim.latmax],...
   [SOCAT_grid.lim.lonmin ...
    SOCAT_grid.lim.lonmax]);
setm(gca,'MapProjection','miller');
set(gca,'fontsize',16);
pcolorm(repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1),...
        repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y),...
        mean(SOCAT_grid.fco2_ave_wtd_detrend,3,'omitnan'));
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(parula(20));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Surface {\itf}_{CO2(SOCAT)} (\muatm)';
cbarrow;
mlabel off;
plabel off;

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/fCO2.png');
close
% clean up
clear
