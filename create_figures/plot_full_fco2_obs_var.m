%% Plot all detrended gridded mean fCO2
% define regions
define_regions_eiwg
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
% title('DJF');
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','south','Position',[0.45 0.2 0.3 0.025]);
colormap(flipud(hot));
caxis([-2.5 82.5]);
c.TickLength = 0;
c.FontWeight = 'bold';
c.FontSize = 10;
c.Label.String = 'Surface {\itf}_{CO2} Variability (\muatm)';
cbarrow;
% plot background
z = std(SOCAT_grid.fco2_ave_wtd_detrend,[],3,'omitnan')';
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,z);
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    z = std(SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,[],3,'omitnan')';
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,z);
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
% figure properties
plot_land('map');
%mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/fCO2_obs_var.png');
close
% clean up
clear
