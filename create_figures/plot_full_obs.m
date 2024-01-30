%% Plot number of all grid cells with data
% define regions
define_regions_eiwg
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','south','Position',[0.45 0.2 0.3 0.025]);
mycolormap = jet(21);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-1 41]);
c.FontSize = 10;
c.FontWeight = 'bold';
c.TickLength = 0;
c.Label.String = 'Total Months Represented';
cbarrow('right');
% plot background
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,SOCAT_grid.num_months');
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,...
        SOCAT_grid.(region{n}).num_months');
    clear SOCAT_grid
end
% plot borders around regions
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% plot land
plot_land('map');
%mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/num_obs.png');
close
% clean up
clear n h r tmp_lon c mycolormap 
