% Evaluate US LME OA indicators against other gridded products

% define LMEs
define_regions_eiwg

%% load SeaFlux
SeaFlux = netcdfreader('evaluation/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc');
seafrac = netcdfreader('evaluation/SeaFlux.v2023.01_seafrac_1982-2022.nc');
% convert longitude
SeaFlux.lon = convert_lon(SeaFlux.lon);
% convert time
time_temp = datenum(1982,1,1+double(SeaFlux.time));
time_temp = datevec(time_temp);
SeaFlux.time = datenum(time_temp(:,1),time_temp(:,2),15);
clear time_temp

%% load CMEMS-LSCE
% CMEMS_LSCE = netcdfreader('/raid/Data/CMEMS_LSCE/CMEMS_LSCE.nc');
% % convert time
% time_temp = datenum(1950,1,1,double(CMEMS_LSCE.time),0,0);
% time_temp = datevec(time_temp);
% CMEMS_LSCE.time = datenum(time_temp(:,1),time_temp(:,2),15);
% % remove unnecessary variables
% vars = fieldnames(CMEMS_LSCE);
% idx = ismember(vars,{'longitude' 'latitude' 'time' 'spco2'});
% CMEMS_LSCE = rmfield(CMEMS_LSCE,vars(~idx));
% clear idx time_temp vars

%% load OceanSODA-ETHZ
% SODA_ETHZ = netcdfreader('/raid/Data/OceanSODA-ETHZ/OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc');
% % convert longitude
% SODA_ETHZ.lon = convert_lon(SODA_ETHZ.lon);
% % convert time
% time_temp = datenum(1982,1,1+double(SODA_ETHZ.time));
% time_temp = datevec(time_temp);
% SODA_ETHZ.time = datenum(time_temp(:,1),time_temp(:,2),15);
% % remove unnecessary variables
% vars = fieldnames(SODA_ETHZ);
% idx = ismember(vars,{'lon' 'lat' 'time' 'spco2'});
% SODA_ETHZ = rmfield(SODA_ETHZ,vars(~idx));

%% load MPI-SOMFFN
% MPI_SOMFFN = netcdfreader('/raid/Data/MPI_SOM_FFN/MPI_SOM_FFN_2022_NCEI_OCADS.nc');
% % convert time
% time_temp = datenum(2000,1,1,0,0,double(MPI_SOMFFN.time));
% time_temp = datevec(time_temp);
% MPI_SOMFFN.time = datenum(time_temp(:,1),time_temp(:,2),15);
% % remove unnecessary variables
% vars = fieldnames(MPI_SOMFFN);
% idx = ismember(vars,{'lon' 'lat' 'time' 'spco2_raw'});
% MPI_SOMFFN = rmfield(MPI_SOMFFN,vars(~idx));

%%
% add pre-allocated pCO2 difference grids
SeaFlux.OAI_grid_mean = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.JENA_MLS_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.MPI_SOMFFN_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.CMEMS_FFNN_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.CSIR_ML6_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.JMA_MLR_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
SeaFlux.NIES_FNN_diff = nan(length(SeaFlux.lon),length(SeaFlux.lat));
datasets = {'JENA_MLS' 'MPI_SOMFFN' 'CMEMS_FFNN' 'CSIR_ML6' 'JMA_MLR' 'NIES_FNN'};
for n = 1:length(region)

    %% remove observations outside general LME limits
    % determine geographic indices
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X)';
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_xmin = find(abs(SeaFlux.lon - min(tmp_lon))==...
        min(abs(SeaFlux.lon - min(tmp_lon))));
    idx_xmax = find(abs(SeaFlux.lon - max(tmp_lon))==...
        min(abs(SeaFlux.lon - max(tmp_lon))));
    idx_ymin = find(abs(SeaFlux.lat - min(tmp_lat))==...
        min(abs(SeaFlux.lat - min(tmp_lat))));
    idx_ymax = find(abs(SeaFlux.lat - max(tmp_lat))==...
        min(abs(SeaFlux.lat - max(tmp_lat))));
    % remove gridded observations outside region
    for d = 1:length(datasets)
       if idx_xmin > idx_xmax
           SeaFlux.(region{n}).(datasets{d}) = ...
               [SeaFlux.(datasets{d})(idx_xmin:end,idx_ymin:idx_ymax,:);...
               SeaFlux.(datasets{d})(1:idx_xmax,idx_ymin:idx_ymax,:)];
       else
           SeaFlux.(region{n}).(datasets{d}) = ...
               SeaFlux.(datasets{d})(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
       end
    end
    % copy other variables
    if idx_xmin > idx_xmax
        SeaFlux.(region{n}).lon = [SeaFlux.lon(idx_xmin:end);SeaFlux.lon(1:idx_xmax)];
    else
        SeaFlux.(region{n}).lon = SeaFlux.lon(idx_xmin:idx_xmax);
    end
    SeaFlux.(region{n}).lat = SeaFlux.lat(idx_ymin:idx_ymax);
    SeaFlux.(region{n}).time = SeaFlux.time;
    % clean up
    clear tmp_lon tmp_lat idx_xmin idx_xmax idx_ymin idx_ymax vars v
   
%     % test plot
%     figure;
%     pcolor(SeaFlux.(region{n}).lon,SeaFlux.(region{n}).lat,...
%         mean(SeaFlux.(region{n}).(datasets{1}),3,'omitnan')');

    %% remove observations outside refined LME limits
    % determine index based on LME
    SeaFlux.(region{n}).idxspc = ...
        nan(length(SeaFlux.(region{n}).lon),length(SeaFlux.(region{n}).lat));
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    SeaFlux.(region{n}).idxspc = ...
        inpolygon(...
        repmat(SeaFlux.(region{n}).lon,1,length(SeaFlux.(region{n}).lat)),...
        repmat(SeaFlux.(region{n}).lat',length(SeaFlux.(region{n}).lon),1),...
        tmp_lon,tmp_lat);
    SeaFlux.(region{n}).idxspc = ...
        repmat(SeaFlux.(region{n}).idxspc,1,1,length(SeaFlux.(region{n}).time));
    % eliminate gridded data outside LME
    for d = 1:length(datasets)
        SeaFlux.(region{n}).(datasets{d})(~SeaFlux.(region{n}).idxspc) = NaN;
    end

%     % test plot
%     figure;
%     pcolor(SeaFlux.(region{n}).lon,SeaFlux.(region{n}).lat,...
%         mean(SeaFlux.(region{n}).(datasets{1}),3,'omitnan')');

    %% load US LME data
    load(['Data/' region{n} '/ML_fCO2.mat']);
    OAI_grid.(region{n}).time = datenum(1998,OAI_grid.(region{n}).month+0.5,15);
    
%     % test plot
%     figure;
%     pcolor(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,...
%         mean(OAI_grid.(region{n}).pCO2,3,'omitnan')'); hold on
%     plot(tmp_lon,tmp_lat)

    %% calculate means of gridded product over time of LME-RFR
    idx_time = SeaFlux.(region{n}).time >= min(OAI_grid.(region{n}).time) & ...
        SeaFlux.(region{n}).time <= max(OAI_grid.(region{n}).time);
    for d = 1:length(datasets)
       SeaFlux.(region{n}).([datasets{d} '_mean']) = ...
           mean(SeaFlux.(region{n}).(datasets{d})(:,:,idx_time),3,'omitnan');
    end

    %% interpolate US LME data onto coarse grid
    OAI_grid.(region{n}).pCO2_mean = mean(OAI_grid.(region{n}).pCO2,3,'omitnan');
    [lon,lat] = meshgrid(SeaFlux.(region{n}).lon,SeaFlux.(region{n}).lat);
    SeaFlux.(region{n}).OAI_grid_mean = ...
        interp2(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).pCO2_mean',lon,lat)';
    for d = 1:length(datasets)
        SeaFlux.(region{n}).([datasets{d} '_diff']) = ...
            SeaFlux.(region{n}).([datasets{d} '_mean'])-...
            SeaFlux.(region{n}).OAI_grid_mean;
    end

%     % test plot
%     figure;
%     pcolor(OAI_grid.(region{n}).(datasets{1}).lon,...
%         OAI_grid.(region{n}).(datasets{1}).lat,...
%         OAI_grid.(region{n}).(datasets{1}).pCO2_mean'); hold on
%     plot_land('xy')
%     plot(tmp_lon,tmp_lat)

    % clean up
    clear OAI_grid

    %% add interpolated values to full SeaFlux grid
    % region-specific indices
    idx_minlat = find(abs(SeaFlux.lat - min(SeaFlux.(region{n}).lat)) == ...
        min(abs(SeaFlux.lat - min(SeaFlux.(region{n}).lat))));
    idx_maxlat = find(abs(SeaFlux.lat - max(SeaFlux.(region{n}).lat)) == ...
        min(abs(SeaFlux.lat - max(SeaFlux.(region{n}).lat))));
    idx_minlon = find(abs(SeaFlux.lon - min(SeaFlux.(region{n}).lon)) == ...
        min(abs(SeaFlux.lon - min(SeaFlux.(region{n}).lon))));
    idx_maxlon = find(abs(SeaFlux.lon - max(SeaFlux.(region{n}).lon)) == ...
        min(abs(SeaFlux.lon - max(SeaFlux.(region{n}).lon))));
    % add LME-RFR mean and differences to SeaFlux grid
    if idx_minlon > idx_maxlon
        SeaFlux.OAI_grid_mean(idx_minlon:end,idx_minlat:idx_maxlat,:) = ...
            SeaFlux.(region{n}).OAI_grid_mean(1:361-idx_minlon,:);
        SeaFlux.OAI_grid_mean(1:idx_maxlon,idx_minlat:idx_maxlat,:) = ...
            SeaFlux.(region{n}).OAI_grid_mean(1:idx_maxlon,:);
        for d = 1:length(datasets)
            SeaFlux.([datasets{d} '_diff'])(idx_minlon:end,idx_minlat:idx_maxlat,:) = ...
                SeaFlux.(region{n}).([datasets{d} '_diff'])(1:361-idx_minlon,:);
            SeaFlux.([datasets{d} '_diff'])(1:idx_maxlon,idx_minlat:idx_maxlat,:) = ...
                SeaFlux.(region{n}).([datasets{d} '_diff'])(1:idx_maxlon,:);
        end
    else
        SeaFlux.OAI_grid_mean(idx_minlon:idx_maxlon,idx_minlat:idx_maxlat,:) = ...
            SeaFlux.(region{n}).OAI_grid_mean;
        for d = 1:length(datasets)
            SeaFlux.([datasets{d} '_diff'])(idx_minlon:idx_maxlon,idx_minlat:idx_maxlat,:) = ...
                SeaFlux.(region{n}).([datasets{d} '_diff']);
        end
    end

end

%% plot all
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
pcolorm(SeaFlux.lat,SeaFlux.lon,mean(SeaFlux.OAI_grid_mean,3,'omitnan')');
plot_land('map');
% figure properties
c=colorbar('location','southoutside');
c.TickLength = 0;
c.Label.String = '';
clim([-50 50]);
colormap(cmocean('balance','pivot',0));
cbarrow;
mlabel off;
% plot differences
pcolorm(SeaFlux.lat,SeaFlux.lon,SeaFlux.mean(JENA_MLS,3,'omitnan')');

contourfm(SeaFlux.lat,SeaFlux.lon,SeaFlux.OAI_grid_mean');
pcolorm(SeaFlux.lat,SeaFlux.lon,SeaFlux.OAI_grid_mean');
pcolorm(SeaFlux.lat,SeaFlux.lon,SeaFlux.JENA_MLS_diff');
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
mlabel off
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/obs.png');
close
% clean up
clear land c mycolormap

%% determine mean and standard deviations of differences

%% Add full LME-RFR fCO2 to full grid
% Load gridded mooring observations
load('Data/socat_gridded','SOCAT_grid');
% Pre-allocate LME RFR grid
US_LME_RFR.lim = SOCAT_grid.lim;
US_LME_RFR.dim = SOCAT_grid.dim;
US_LME_RFR.lon = SOCAT_grid.lon;
US_LME_RFR.lat = SOCAT_grid.lat;
US_LME_RFR.month = SOCAT_grid.month;
US_LME_RFR.year = SOCAT_grid.year;
US_LME_RFR.(datasets{d})
clear SOCAT_grid m
% define LMEs
define_regions_eiwg
% 
for n = 1:length(region)
    % load gridded fCO2
    load(['Data/' region{n} '/ML_fCO2_gridded_comparison'],'OAI_grid');
    % Add fCO2 to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);

    US_LME_RFR.fCO2(idx_lon,idx_lat,:) = OAI_grid.(region{n}).fCO2;
    clear idx_lon idx_lat OAI_grid n
end


%% Pre-allocate LME RFR grid
US_LME_RFR.lim = SOCAT_grid.lim;
US_LME_RFR.dim = SOCAT_grid.dim;
US_LME_RFR.lon = SOCAT_grid.lon;
US_LME_RFR.lat = SOCAT_grid.lat;
US_LME_RFR.month = SOCAT_grid.month;
US_LME_RFR.year = SOCAT_grid.year;
US_LME_RFR.fCO2_no_moorings = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
US_LME_RFR.fCO2_full = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
clear SOCAT_grid m
