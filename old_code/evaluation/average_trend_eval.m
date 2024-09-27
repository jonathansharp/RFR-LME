% Evaluate US LME OA indicators against gridded products

%% load US LME data
date = '02-Aug-2023';
LME_RFR = netcdfreader(['Data/US_LME_RFR_Inds_' date '.nc']);
LME_RFR_u = netcdfreader(['Data/US_LME_RFR_Inds_Uncer_' date '.nc']);

%% load CMEMS-LSCE
CMEMS_LSCE = netcdfreader('/raid/Data/CMEMS_LSCE/CMEMS_LSCE.nc');
CMEMS_LSCE.time = datenum(1950,1,1,double(CMEMS_LSCE.time),0,0);

%% load OceanSODA-ETHZ
SODA_ETHZ = netcdfreader('/raid/Data/OceanSODA-ETHZ/OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc');
SODA_ETHZ.time = datenum(1982,1,1+double(SODA_ETHZ.time));
SODA_ETHZ.lon = convert_lon(SODA_ETHZ.lon);

%% load MPI-SOMFFN
MPI_SOMFFN = netcdfreader('/raid/Data/MPI_SOM_FFN/MPI_SOM_FFN_2022_NCEI_OCADS.nc');
MPI_SOMFFN.time = datenum(2000,1,1,0,0,double(MPI_SOMFFN.time));

%% plot regional average trends
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
clim([0.5 2.5]);
colormap(cmocean('amp',18));
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itf}CO_{2} Trend (\muatm yr^{-1})';
% plot land
plot_land('map');
mlabel off
% plot regions
for n = 1:length(region)
    % load
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % calculate area-weighted time series
    OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
    area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
    for t = 1:OAI_grid.(region{n}).dim.z
        % remove ice-filled cells
        area_weights(isnan(OAI_grid.(region{n}).fCO2(:,:,t))) = NaN;
        OAI_grid.(region{n}).fCO2_dom_mean(t) = ...
            squeeze(sum(sum(OAI_grid.(region{n}).fCO2(:,:,t).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    end
    % calculate trend
    [yf,yr,x] = ...
        leastsq2(OAI_grid.(region{n}).month,...
        OAI_grid.(region{n}).fCO2_dom_mean,0,2,[6 12]);
    % replicate trend to grid
    trend_grid = nan(length(OAI_grid.(region{n}).lon),length(OAI_grid.(region{n}).lat));
    trend_grid(OAI_grid.(region{n}).idxspc(:,:,1)) = x(2)*12;
    % plot
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        trend_grid',0.5:0.125:2.5,'LineStyle','none');
    clear vars_grid z
    % 
    disp(['RFR-LME Trend (' region{n} ') = ' ...
        num2str(round(mean(x(2)*12,3,'omitnan'),1)) ' uatm/yr']);
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
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Average_trend_map.png');

