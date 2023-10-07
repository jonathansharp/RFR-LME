% Evaluate US LME OA indicators against other gridded products

%% define LMEs
define_regions_eiwg

%% define datasets
datasets = {'JENA_MLS' 'MPI_SOMFFN' 'CMEMS_FFNN' 'CSIR_ML6' 'JMA_MLR' 'NIES_FNN'};

%% load US LME data
date = '02-Aug-2023';
LME_RFR = netcdfreader(['Data/US_LME_RFR_Inds_' date '.nc']);
clear date

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
LME_RFR.pCO2_mean = mean(LME_RFR.pCO2,3,'omitnan');

%% calculate means of gridded products over time and fill with Landschutzer climatology
% index to timespan of RFR-LME
idx_time = SeaFlux.time >= min(LME_RFR.Time) & SeaFlux.time <= max(LME_RFR.Time);
% determine long-term mean for each dataset
for d = 1:length(datasets)
   SeaFlux.([datasets{d} '_mean']) = ...
       mean(SeaFlux.(datasets{d})(:,:,idx_time),3,'omitnan');
end
% fill empty cells in each dataset with climatology
SeaFluxFiller = netcdfreader('evaluation/SeaFlux_v2021.04_spco2_filler_1990-2019.nc');
mean_filler = mean(SeaFluxFiller.spco2_scaled_climatology,3,'omitnan');
for d = 1:length(datasets)
    idx = isnan(SeaFlux.([datasets{d} '_mean']));
    SeaFlux.([datasets{d} '_mean'])(idx) = mean_filler(idx);
end
clear idx_time SeaFluxFiller mean_filler idx d

%% interpolate RFR-LME to coarse grid and calculate differences between gridded products and RFR-LME
[lon,lat] = meshgrid(SeaFlux.lon,SeaFlux.lat);
SeaFlux.RFR_LME_mean = interp2(LME_RFR.Lon,LME_RFR.Lat,LME_RFR.pCO2_mean',lon,lat)';
for d = 1:length(datasets)
    SeaFlux.([datasets{d} '_diff']) = SeaFlux.([datasets{d} '_mean'])-SeaFlux.RFR_LME_mean;
end
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
    clim([295 475]);
    colormap(parula);
    cbarrow;
    mlabel off;
    % plot differences
    pcolorm(SeaFlux.lat-0.5,[SeaFlux.lon;SeaFlux.lon(end)+1]-0.5,...
        [SeaFlux.([datasets{d} '_mean']);SeaFlux.([datasets{d} '_mean'])(end,:)]');
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
    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/gridded_diffs_' datasets{d} '.png']);
    close
    clear c n tmp_lon tmp_lat
end
clear d

%% Determine area of each grid cell
seafrac = netcdfreader('evaluation/SeaFlux.v2023.01_seafrac_1982-2022.nc');
area_km2_temp = ...
(((repmat(SeaFlux.lat',length(SeaFlux.lon),1) + 0.5) - ...
    (repmat(SeaFlux.lat',length(SeaFlux.lon),1) - 0.5)) .* 110.574) .* ... % latitude distance
(((repmat(SeaFlux.lon,1,length(SeaFlux.lat)) + 0.5) - ...
    (repmat(SeaFlux.lon,1,length(SeaFlux.lat)) - 0.5)) .* ...
    111.320.*cosd(repmat(SeaFlux.lat',length(SeaFlux.lon),1))); % longitude distance
SeaFlux.area_km2 = area_km2_temp.*seafrac.seafrac;
clear area_km2_temp seafrac

%% calculate are-weighted means and standard deviations
for d = 1:length(datasets)
    idx = ~isnan(SeaFlux.([datasets{d} '_diff']));
    weighted_mean = sum(SeaFlux.([datasets{d} '_diff'])(idx).*...
        SeaFlux.area_km2(idx))./sum(SeaFlux.area_km2(idx))
    weighted_std = ...
        std(SeaFlux.([datasets{d} '_diff'])(idx),SeaFlux.area_km2(idx))
end

%% clear up
clear