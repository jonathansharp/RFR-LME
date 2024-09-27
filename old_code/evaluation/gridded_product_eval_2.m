% Evaluate US LME OA indicators against other gridded products

%% define LMEs
define_regions_eiwg

%% define datasets
datasets = {'JENA_MLS' 'MPI_SOMFFN' 'CMEMS_FFNN' 'CSIR_ML6' 'JMA_MLR' 'NIES_FNN'};

%% load US LME data

% RFR-LME file data!
date = '08-May-2024';
%date = '06-Oct-2023';
new = 1;

% load data
if new == 1
% New NetCDF files
RFR_LME.Lat = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'lat');
RFR_LME.Lon = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'lon');
RFR_LME.Time = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'time');
RFR_LME.pCO2 = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'pco2');
RFR_LME.fCO2 = ncread(['Data/NetCDFs_' date '/US_RFR_LME_fCO2.nc'],'fco2');
RFR_LME.TA = ncread(['Data/NetCDFs_' date '/US_RFR_LME_TA.nc'],'ta');
RFR_LME.DIC = ncread(['Data/NetCDFs_' date '/US_RFR_LME_DIC.nc'],'dic');
RFR_LME.pH = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pH.nc'],'ph');
RFR_LME.OmA = ncread(['Data/NetCDFs_' date '/US_RFR_LME_OmA.nc'],'om_a');
RFR_LME.OmC = ncread(['Data/NetCDFs_' date '/US_RFR_LME_OmC.nc'],'om_c');
RFR_LME.H = ncread(['Data/NetCDFs_' date '/US_RFR_LME_H.nc'],'h');
RFR_LME.CO3 = ncread(['Data/NetCDFs_' date '/US_RFR_LME_CO3.nc'],'co3');
RFR_LME.RF = ncread(['Data/NetCDFs_' date '/US_RFR_LME_RF.nc'],'rf');
else
% Old NetCDF files
RFR_LME.Lat = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Lat');
RFR_LME.Lon = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Lon');
RFR_LME.Time = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Time');
RFR_LME.pCO2 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'pCO2');
RFR_LME.fCO2 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'fCO2');
RFR_LME.TA = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'TA');
RFR_LME.DIC = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'DIC');
RFR_LME.pH = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'pH');
RFR_LME.OmA = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'OmA');
RFR_LME.OmC = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'OmC');
RFR_LME.H = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'H');
RFR_LME.CO3 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'CO3');
RFR_LME.RF = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'RF');
end


%% load SeaFlux
SeaFlux = netcdfreader('evaluation/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc');
% convert longitude
SeaFlux.lon = convert_lon(SeaFlux.lon);
% convert time
time_temp = datenum(1982,1,1+double(SeaFlux.time));
time_temp = datevec(time_temp);
SeaFlux.time = datenum(time_temp(:,1),time_temp(:,2),15);
clear time_temp

%% calculate mean of RFR-LME  over time
idx_time = RFR_LME.Time >= min(SeaFlux.time) & RFR_LME.Time <= max(SeaFlux.time);
RFR_LME.pCO2_mean = mean(RFR_LME.pCO2(:,:,idx_time),3,'omitnan');

%% calculate means of gridded products over time and fill with Landschutzer climatology
% index to timespan of RFR-LME
idx_time = SeaFlux.time >= min(RFR_LME.Time) & SeaFlux.time <= max(RFR_LME.Time);
% determine long-term mean for each dataset
for d = 1:length(datasets)
   SeaFlux.([datasets{d} '_mean']) = ...
       mean(SeaFlux.(datasets{d})(:,:,idx_time),3,'omitnan');
   SeaFlux.([datasets{d} '_mean_filled']) = ...
       mean(SeaFlux.(datasets{d})(:,:,idx_time),3,'omitnan');
end
% fill empty cells in each dataset with climatology
SeaFluxFiller = netcdfreader('evaluation/SeaFlux_v2021.04_spco2_filler_1990-2019.nc');
idx_time = datenum(1982,1,15) + double(SeaFluxFiller.time) >= min(RFR_LME.Time) & ...
    datenum(1982,1,15) + double(SeaFluxFiller.time) <= max(RFR_LME.Time);
mean_filler = mean(SeaFluxFiller.spco2_scaled_climatology(:,:,idx_time),3,'omitnan');
for d = 1:length(datasets)
    idx = isnan(SeaFlux.([datasets{d} '_mean']));
    SeaFlux.([datasets{d} '_mean_filled'])(idx) = mean_filler(idx);
end
clear idx_time SeaFluxFiller mean_filler idx d
SeaFlux.ensemble_mean = mean(cat(3,SeaFlux.JENA_MLS_mean,...
    SeaFlux.MPI_SOMFFN_mean,SeaFlux.CMEMS_FFNN_mean,...
    SeaFlux.CSIR_ML6_mean,SeaFlux.JMA_MLR_mean,...
    SeaFlux.NIES_FNN_mean),3,'omitnan');
SeaFlux.ensemble_mean_filled = mean(cat(3,SeaFlux.JENA_MLS_mean_filled,...
    SeaFlux.MPI_SOMFFN_mean_filled,SeaFlux.CMEMS_FFNN_mean_filled,...
    SeaFlux.CSIR_ML6_mean_filled,SeaFlux.JMA_MLR_mean_filled,...
    SeaFlux.NIES_FNN_mean_filled),3,'omitnan');

%% interpolate RFR-LME to coarse grid and calculate differences between gridded products and RFR-LME
[lon,lat] = meshgrid(SeaFlux.lon,SeaFlux.lat);
SeaFlux.RFR_LME_mean = interp2(RFR_LME.Lon,RFR_LME.Lat,RFR_LME.pCO2_mean',lon,lat)';
for d = 1:length(datasets)
    SeaFlux.([datasets{d} '_diff']) = SeaFlux.([datasets{d} '_mean'])-SeaFlux.RFR_LME_mean;
    SeaFlux.([datasets{d} '_diff_filled']) = SeaFlux.([datasets{d} '_mean_filled'])-SeaFlux.RFR_LME_mean;
end
% compute difference with SeaFlux ensemble mean
SeaFlux.ensemble_diff = SeaFlux.ensemble_mean-SeaFlux.RFR_LME_mean;
SeaFlux.ensemble_diff_filled = SeaFlux.ensemble_mean_filled-SeaFlux.RFR_LME_mean;
datasets = [datasets {'ensemble'}];
clear lon lat d

%% plot gridded products within each LME
for d = 1:length(datasets)
    % initialize figure
    figure('visible','on'); box on; hold on;
    worldmap([-18 82],[140 302]);
    setm(gca,'MapProjection','robinson','MLabelParallel','south');
    set(gcf,'position',[100 100 900 600]);
    set(gca,'fontsize',16);
    % figure properties
    c=colorbar('location','southoutside');
    c.TickLength = 0;
    c.Label.String = ['Sea Surface {\itp}CO_{2(' strrep(datasets{d},'_','-') ')} (\muatm)'];
    clim([295 500]);
    colormap(flipud(slanCM('romao')));
    cbarrow;
    mlabel off;
    % plot differences
    pcolorm(SeaFlux.lat-0.5,[SeaFlux.lon;SeaFlux.lon(end)+1]-0.5,...
        [SeaFlux.([datasets{d} '_mean']);SeaFlux.([datasets{d} '_mean'])(end,:)]');
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    plot_land('map');
    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/gridded_mean_' datasets{d} '.png']);
    close
    clear c n tmp_lon tmp_lat
end
clear d

%% plot differences between gridded products and RFR-LME
for d = 1:length(datasets)
    % initialize figure
    figure('visible','on'); box on; hold on;
    worldmap([-18 82],[140 302]);
    setm(gca,'MapProjection','robinson','MLabelParallel','south');
    set(gcf,'position',[100 100 900 600]);
    set(gca,'fontsize',16);
    % figure properties
    c=colorbar('location','southoutside');
    c.TickLength = 0;
    c.Label.String = ['\Delta{\itp}CO_{2(' strrep(datasets{d},'_','-') ')} (\muatm)'];
    clim([-50 50]);
    colormap(cmocean('balance','pivot',0));
    cbarrow;
    mlabel off;
    % plot differences
    pcolorm(SeaFlux.lat-0.5,[SeaFlux.lon;SeaFlux.lon(end)+1]-0.5,...
        [SeaFlux.([datasets{d} '_diff']);SeaFlux.([datasets{d} '_diff'])(end,:)]');
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    plot_land('map');
    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/gridded_diffs_' datasets{d} '.png']);
    close
    clear c n tmp_lon tmp_lat
end
clear d

%% Determine area of each grid cell
seafrac = netcdfreader('evaluation/SeaFlux.v2023.02_seafrac_1982-2022.nc');
area_km2_temp = ...
(((repmat(SeaFlux.lat',length(SeaFlux.lon),1) + 0.5) - ...
    (repmat(SeaFlux.lat',length(SeaFlux.lon),1) - 0.5)) .* 110.574) .* ... % latitude distance
(((repmat(SeaFlux.lon,1,length(SeaFlux.lat)) + 0.5) - ...
    (repmat(SeaFlux.lon,1,length(SeaFlux.lat)) - 0.5)) .* ...
    111.320.*cosd(repmat(SeaFlux.lat',length(SeaFlux.lon),1))); % longitude distance
SeaFlux.area_km2 = area_km2_temp.*seafrac.seafrac;
clear area_km2_temp seafrac

%% calculate area-weighted means and standard deviations of unfilled products
for d = 1:length(datasets)
    idx = ~isnan(SeaFlux.([datasets{d} '_diff']));
    weighted_mean = sum(SeaFlux.([datasets{d} '_diff'])(idx).*...
        SeaFlux.area_km2(idx))./sum(SeaFlux.area_km2(idx));
    weighted_std = std(SeaFlux.([datasets{d} '_diff'])(idx),SeaFlux.area_km2(idx));
    disp([datasets{d} ' (unfilled) = ' num2str(round(weighted_mean,1)) ...
        ' +/- ' num2str(round(weighted_std,1))]);
end

%% calculate area-weighted means and standard deviations of filled products
for d = 1:length(datasets)
    idx = ~isnan(SeaFlux.([datasets{d} '_diff_filled']));
    weighted_mean = sum(SeaFlux.([datasets{d} '_diff_filled'])(idx).*...
        SeaFlux.area_km2(idx))./sum(SeaFlux.area_km2(idx));
    weighted_std = std(SeaFlux.([datasets{d} '_diff_filled'])(idx),SeaFlux.area_km2(idx));
    disp([datasets{d} ' (filled) = ' num2str(round(weighted_mean,1)) ...
        ' +/- ' num2str(round(weighted_std,1))]);
end

%% clear up
clear