%% Plot number of climatological grid cells with data
% load SOCAT grid
load('Data/socat_gridded_no_moorings','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
mycolormap = jet(13);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-0.5 12.5]);
c.TickLength = 0;
c.Label.String = 'Months of Year Represented';
% plot background
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,SOCAT_grid.num_months_clim');
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,...
        SOCAT_grid.(region{n}).num_months_clim');
    clear SOCAT_grid
end
% plot borders around regions
plot_lme_borders(region,lme_shape,lme_idx);
% plot land
plot_land('map');
mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/num_obs_clim.png');
close
% clean up
clear n h r tmp_lon c mycolormap 
