% Compare LME-RFR to mooring observations
% 
% This script compares pCO2 predictions in US LMEs to gridded pCO2 at
% mooring observation sites, obtained from SOCAT
% 
% Written by J.D. Sharp: 8/1/23
% Last updated by J.D. Sharp: 10/18/24

%% Load gridded mooring observations
load('Data/socat_gridded_2023_moorings_only','SOCAT_grid');
% define mooring indices by grid cells that include data
[a,b] = find(~isnan(mean(SOCAT_grid.fco2_ave_wtd,3,'omitnan')));

%% extract mooring timeseries (only when mooring has >36 months represented)
cnt=1; % set counter to one
for m = 1:length(a)
    mooring.(['m' num2str(m)]).fCO2 = ...
        squeeze(SOCAT_grid.fco2_ave_wtd(a(m),b(m),:));
    mooring.(['m' num2str(m)]).SST = ...
        squeeze(SOCAT_grid.sst_ave_wtd(a(m),b(m),:));
    mooring.(['m' num2str(m)]).SSS = ...
        squeeze(SOCAT_grid.sss_ave_wtd(a(m),b(m),:));
    C=CO2SYS(2400,mooring.(['m' num2str(m)]).fCO2,1,5,mooring.(['m' num2str(m)]).SSS,...
        mooring.(['m' num2str(m)]).SST,NaN,0,NaN,0,0,0,0,1,10,1,2,2);
    pCO2 = C(:,4); pCO2(pCO2==-999)=NaN;
    mooring.(['m' num2str(m)]).pCO2 = pCO2;
    mooring.(['m' num2str(m)]).lon = SOCAT_grid.lon(a(m));
    mooring.(['m' num2str(m)]).lat = SOCAT_grid.lat(b(m));
    mooring.(['m' num2str(m)]).time = datenum(1998,SOCAT_grid.month+0.5,15);
    % average together the two WHOTS grid cells
    if mooring.(['m' num2str(m)]).lon == 202.125 && ...
            mooring.(['m' num2str(m)]).lat > 22.5 && ...
            mooring.(['m' num2str(m)]).lat < 23
        temp(:,cnt) = mooring.(['m' num2str(m)]).pCO2;
        WHOTSm(cnt) = m;
        cnt=cnt+1;
    end
    % remove moorings that aren't in LMEs or don't have enough obs
    if sum(~isnan(mooring.(['m' num2str(m)]).pCO2)) < 36
        mooring = rmfield(mooring,(['m' num2str(m)]));
    end
end
mooring.(['m' num2str(WHOTSm(1))]).pCO2 = mean(temp,2,'omitnan');
mooring = rmfield(mooring,(['m' num2str(WHOTSm(2))]));
moor_nums = fieldnames(mooring);
% clean up
clear a b m cnt temp WHOTSm

%% Pre-allocate LME RFR grid
US_LME_RFR.lim = SOCAT_grid.lim;
US_LME_RFR.dim = SOCAT_grid.dim;
US_LME_RFR.lon = SOCAT_grid.lon;
US_LME_RFR.lat = SOCAT_grid.lat;
US_LME_RFR.month = SOCAT_grid.month;
US_LME_RFR.year = SOCAT_grid.year;
US_LME_RFR.pCO2_no_moorings = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
US_LME_RFR.pCO2_full = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
clear SOCAT_grid m

%% Add mooring-excluded LME-RFR pCO2 to full grid
define_regions_eiwg
for n = 1:length(region)
    % load gridded pCO2
    load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');
    % Add pCO2 to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);
    US_LME_RFR.pCO2_no_moorings(idx_lon,idx_lat,:) = OAI_grid.(region{n}).pCO2;
%     US_LME_RFR.upco2(idx_lon,idx_lat,:) = OAI_grid.(region{m}).upCO2;
    % clean up
    clear idx_lon idx_lat OAI_grid n
end
%% Add full LME-RFR pCO2 to full grid
define_regions_eiwg
for n = 1:length(region)
    % load gridded fCO2
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % Add pCO2 to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);
    US_LME_RFR.pCO2_full(idx_lon,idx_lat,:) = OAI_grid.(region{n}).pCO2;
%     US_LME_RFR.upco2(idx_lon,idx_lat,:) = OAI_grid.(region{m}).upCO2;
    % clean up
    clear idx_lon idx_lat OAI_grid n
end

%% extract time series from LME-RFR data
for m = 1:length(moor_nums)
    londiffs = sort(abs(US_LME_RFR.lon - mooring.(moor_nums{m}).lon));
    latdiffs = sort(abs(US_LME_RFR.lat - mooring.(moor_nums{m}).lat));
    lon_idx = find(abs(US_LME_RFR.lon - mooring.(moor_nums{m}).lon) == londiffs(1));
    lat_idx = find(abs(US_LME_RFR.lat - mooring.(moor_nums{m}).lat) == latdiffs(1));
    % try moving lat a bit
    if isnan(US_LME_RFR.pCO2_full(lon_idx,lat_idx,1))
        lat_idx = find(abs(US_LME_RFR.lat - mooring.(moor_nums{m}).lat) == latdiffs(2));
    end
    % move lat back, then try moving lon a bit
    if isnan(US_LME_RFR.pCO2_full(lon_idx,lat_idx,1))
        lat_idx = find(abs(US_LME_RFR.lat - mooring.(moor_nums{m}).lat) == latdiffs(1));
        lon_idx = find(abs(US_LME_RFR.lon - mooring.(moor_nums{m}).lon) == londiffs(2));
    end
    mooring.(moor_nums{m}).US_LME_RFR_pCO2_full = ...
        squeeze(US_LME_RFR.pCO2_full(lon_idx(1),lat_idx(1),:));
    mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings = ...
        squeeze(US_LME_RFR.pCO2_no_moorings(lon_idx(1),lat_idx(1),:));
    mooring.(moor_nums{m}).US_LME_RFR_datetime = ...
        datenum(1998,US_LME_RFR.month+0.5,15);
    mooring.(moor_nums{m}).US_LME_RFR_lat = US_LME_RFR.lat(lat_idx(1));
    mooring.(moor_nums{m}).US_LME_RFR_lon = US_LME_RFR.lon(lon_idx(1));
end

%% load SeaFlux
SeaFlux = netcdfreader('evaluation/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc');
fco2atm = netcdfreader('evaluation/SeaFlux.v2023.02_fco2atm_1982-2022.nc');
seafrac = netcdfreader('evaluation/SeaFlux.v2023.02_seafrac_1982-2022.nc');
idx_time = fco2atm.time <= max(SeaFlux.time);
fco2atm.fco2atm = fco2atm.fco2atm(:,:,idx_time);
% convert time
SeaFlux.time = datenum(1982,1,1+double(SeaFlux.time));
% convert longitude
SeaFlux.lon = convert_lon(SeaFlux.lon);

%% load CMEMS-LSCE
CMEMS_LSCE = netcdfreader('/raid/Data/CMEMS_LSCE/CMEMS_LSCE.nc');
CMEMS_LSCE.time = datenum(1950,1,1,double(CMEMS_LSCE.time),0,0);
vars = fieldnames(CMEMS_LSCE);
idx = ismember(vars,{'longitude' 'latitude' 'time' 'spco2'});
CMEMS_LSCE = rmfield(CMEMS_LSCE,vars(~idx));

%% load OceanSODA-ETHZ
SODA_ETHZ = netcdfreader('/raid/Data/OceanSODA-ETHZ/OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc');
SODA_ETHZ.time = datenum(1982,1,1+double(SODA_ETHZ.time));
SODA_ETHZ.lon = convert_lon(SODA_ETHZ.lon);
vars = fieldnames(SODA_ETHZ);
idx = ismember(vars,{'lon' 'lat' 'time' 'spco2'});
SODA_ETHZ = rmfield(SODA_ETHZ,vars(~idx));

%% load MPI-SOMFFN
MPI_SOMFFN = netcdfreader('/raid/Data/MPI_SOM_FFN/MPI_SOM_FFN_2022_NCEI_OCADS.nc');
MPI_SOMFFN.time = datenum(2000,1,1,0,0,double(MPI_SOMFFN.time));
vars = fieldnames(MPI_SOMFFN);
idx = ismember(vars,{'lon' 'lat' 'time' 'spco2_raw'});
MPI_SOMFFN = rmfield(MPI_SOMFFN,vars(~idx));

%% extract time series from SeaFlux data
for m = 1:length(moor_nums)
    % get initial lat and lon indices
    londiffs = sort(abs(SeaFlux.lon - mooring.(moor_nums{m}).lon));
    latdiffs = sort(abs(SeaFlux.lat - mooring.(moor_nums{m}).lat));
    lon_idx_init = find(abs(SeaFlux.lon - mooring.(moor_nums{m}).lon) == londiffs(1));
    lat_idx_init = find(abs(SeaFlux.lat - mooring.(moor_nums{m}).lat) == latdiffs(1));
    %% JENA_MLS:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.JENA_MLS(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).JENA_MLS_pCO2 = ...
        squeeze(SeaFlux.JENA_MLS(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).JENA_MLS_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);

    %% CSIR_ML6:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.CSIR_ML6(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).CSIR_ML6_pCO2 = ...
        squeeze(SeaFlux.CSIR_ML6(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).CSIR_ML6_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);
    %% JMA_MLR:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.JMA_MLR(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).JMA_MLR_pCO2 = ...
        squeeze(SeaFlux.JMA_MLR(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).JMA_MLR_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);
    %% NIES_FNN:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.NIES_FNN(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).NIES_FNN_pCO2 = ...
        squeeze(SeaFlux.NIES_FNN(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).NIES_FNN_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);
    %% CMEMS-LSCE:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.CMEMS_FFNN(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).CMEMS_LSCE_pCO2 = ...
        squeeze(SeaFlux.CMEMS_FFNN(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).CMEMS_LSCE_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);
    %% MPI-SOMFFN:
    % adjust coordinates until nearest grid cell with a value is found
    [lons,lats] = ndgrid(SeaFlux.lon,SeaFlux.lat);
    lons = lons(:); lats = lats(:);
    lon_idx = lon_idx_init; lat_idx = lat_idx_init;
    while isnan(mean(SeaFlux.MPI_SOMFFN(lon_idx,lat_idx,:),'omitnan'))
       idx = dsearchn([lons(:) lats(:)],...
            [mooring.(moor_nums{m}).lon mooring.(moor_nums{m}).lat]);
       lon_idx = find(SeaFlux.lon == lons(idx));
       lat_idx = find(SeaFlux.lat == lats(idx));
       lons(idx) = []; lats(idx) = [];
    end
    mooring.(moor_nums{m}).MPI_SOMFFN_pCO2 = ...
        squeeze(SeaFlux.MPI_SOMFFN(lon_idx,lat_idx,:));
    datetime_temp = datevec(SeaFlux.time);
    mooring.(moor_nums{m}).MPI_SOMFFN_datetime = ...
        datenum(datetime_temp(:,1),datetime_temp(:,2),15);
end

%% extract time series from CMEMS-LSCE data
% for m = 1:length(moor_nums)
%     londiffs = sort(abs(CMEMS_LSCE.longitude - mooring.(moor_nums{m}).lon));
%     latdiffs = sort(abs(CMEMS_LSCE.latitude - mooring.(moor_nums{m}).lat));
%     lon_idx = find(abs(CMEMS_LSCE.longitude - mooring.(moor_nums{m}).lon) == londiffs(1));
%     lat_idx = find(abs(CMEMS_LSCE.latitude - mooring.(moor_nums{m}).lat) == latdiffs(1));
%     % try moving lat a bit
%     if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1))
%         lat_idx = find(abs(CMEMS_LSCE.latitude - mooring.(moor_nums{m}).lat) == latdiffs(2));
%     end
%     % try moving lon a bit
%     if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1))
%         lon_idx = find(abs(CMEMS_LSCE.longitude - mooring.(moor_nums{m}).lon) == londiffs(2));
%     end
%     mooring.(moor_nums{m}).CMEMS_LSCE_pCO2 = ...
%         squeeze(CMEMS_LSCE.spco2(lon_idx,lat_idx,:));
%     datetime_temp = datevec(CMEMS_LSCE.time);
%     mooring.(moor_nums{m}).CMEMS_LSCE_datetime = ...
%         datenum(datetime_temp(:,1),datetime_temp(:,2),15);
% end

%% extract time series from OceanSODA-ETHZ data
% for m = 1:length(moor_nums)
%     londiffs = sort(abs(SODA_ETHZ.lon - mooring.(moor_nums{m}).lon));
%     latdiffs = sort(abs(SODA_ETHZ.lat - mooring.(moor_nums{m}).lat));
%     lon_idx = find(abs(SODA_ETHZ.lon - mooring.(moor_nums{m}).lon) == londiffs(1));
%     lat_idx = find(abs(SODA_ETHZ.lat - mooring.(moor_nums{m}).lat) == latdiffs(1));
%     % try moving lat a bit
%     if isnan(SODA_ETHZ.spco2(lon_idx,lat_idx,1))
%         lat_idx = find(abs(SODA_ETHZ.lat - mooring.(moor_nums{m}).lat) == latdiffs(2));
%     end
%     % try moving lon a bit
%     if isnan(SODA_ETHZ.spco2(lon_idx,lat_idx,1))
%         lon_idx = find(abs(SODA_ETHZ.lon - mooring.(moor_nums{m}).lon) == londiffs(2));
%     end
%     mooring.(moor_nums{m}).SODA_ETHZ_pCO2 = ...
%         squeeze(SODA_ETHZ.spco2(lon_idx,lat_idx,:));
%     datetime_temp = datevec(SODA_ETHZ.time);
%     mooring.(moor_nums{m}).SODA_ETHZ_datetime = ...
%         datenum(datetime_temp(:,1),datetime_temp(:,2),15);
% end

%% extract time series from MPI-SOMFFN data
% for m = 1:length(moor_nums)
%     londiffs = sort(abs(MPI_SOMFFN.lon - mooring.(moor_nums{m}).lon));
%     latdiffs = sort(abs(MPI_SOMFFN.lat - mooring.(moor_nums{m}).lat));
%     lon_idx = find(abs(MPI_SOMFFN.lon - mooring.(moor_nums{m}).lon) == londiffs(1));
%     lat_idx = find(abs(MPI_SOMFFN.lat - mooring.(moor_nums{m}).lat) == latdiffs(1));
%     % try moving lat a bit
%     if isnan(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,1))
%         lat_idx = find(abs(MPI_SOMFFN.lat - mooring.(moor_nums{m}).lat) == latdiffs(2));
%     end
%     % try moving lon a bit
%     if isnan(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,1))
%         lon_idx = find(abs(MPI_SOMFFN.lon - mooring.(moor_nums{m}).lon) == londiffs(2));
%     end
%     mooring.(moor_nums{m}).MPI_SOMFFN_pCO2 = ...
%         squeeze(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,:));
%     datetime_temp = datevec(MPI_SOMFFN.time);
%     mooring.(moor_nums{m}).MPI_SOMFFN_datetime = ...
%         datenum(datetime_temp(:,1),datetime_temp(:,2),15);
% end


%% plot mooring timeseries comparisons
clrs = cbrewer('qual','Set1',8);
datasets = {'US_LME_RFR' 'US_LME_RFR_NM' 'CMEMS_LSCE' 'NIES_FNN' ...
    'MPI_SOMFFN' 'JENA_MLS' 'CSIR_ML6' 'JMA_MLR'};
for m = 1:length(moor_nums)
    % create plot
    f=figure('visible','off'); hold on;
    set(gca,'fontsize',14);
    f.Position(3) = 2*f.Position(3);
    title([num2str(mooring.(moor_nums{m}).US_LME_RFR_lat) char(176) 'N, ' ...
        num2str(mooring.(moor_nums{m}).US_LME_RFR_lon) char(176) 'E'],'fontweight','bold');
    p.p1=plot(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).US_LME_RFR_pCO2_full,...
        '-','color',clrs(1,:),'linewidth',2);
    p.p2=plot(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings,...
        '--','color',clrs(2,:),'linewidth',2);
    for d = 3:8
        p.(['p' num2str(d)])=...
            plot(mooring.(moor_nums{m}).([datasets{d} '_datetime']),...
            mooring.(moor_nums{m}).([datasets{d} '_pCO2']),...
            '-','color',[clrs(d,:) 0.4],'linewidth',1);
    end
    s1=scatter(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).pCO2,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    % add figure properties
    xlim([min(mooring.(moor_nums{m}).US_LME_RFR_datetime) ...
        max(mooring.(moor_nums{m}).US_LME_RFR_datetime)]);
    ylabel('{\itp}CO_{2} (\muatm)','fontsize',16);
    datetick('x','yyyy','keeplimits');
    if mooring.(moor_nums{m}).lat == 59.875
        ylim([150 500]);
    end
    legend([p.p1 p.p2 p.p3 p.p4 p.p5 p.p6 p.p7 p.p8],...
        {'LME-RFR' 'LME-RFR_{NM}' 'CMEMS-LSCE' 'NIES-FNN' ...
        'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
        'Location','northwest','NumColumns',4);
    % export figure
    % exportgraphics(gcf,['Figures/moorings/m' moor_nums{m} '.png']);
    export_fig(gcf,['Figures/moorings/m' moor_nums{m} '.png'],'-transparent');
    close
end

%% calculate and save error and summary stats
datasets = {'US_LME_RFR' 'US_LME_RFR_NM' 'CMEMS_LSCE' 'NIES_FNN' ...
    'MPI_SOMFFN' 'JENA_MLS' 'CSIR_ML6' 'JMA_MLR'};
% pre-allocate stats
mean_err = nan(7,length(moor_nums));
std_err = nan(7,length(moor_nums));
med_err = nan(7,length(moor_nums));
rmse_err = nan(7,length(moor_nums));
iqr_err = nan(7,length(moor_nums));
resid_corr = nan(7,length(moor_nums));
amp_diff = nan(7,length(moor_nums));
% loop through moorings
for m = 1:length(moor_nums)
    % fit mooring observations to timeseries
    [yf,mooring.(moor_nums{m}).pCO2_resids,x] = ...
        leastsq2(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).pCO2,datenum(1998,1,15),...
        2,[365.25/2 365.25]);
    moor_amp = sqrt(x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2);
    % loop through datasets
    for d = 1:8
        % log dates of gridded values (with month and year)
        if d == 1 || d == 2
            grid_date = datevec(mooring.(moor_nums{m}).US_LME_RFR_datetime);
        else
            grid_date = datevec(mooring.(moor_nums{m}).([datasets{d} '_datetime']));
        end
        grid_date = datenum(grid_date(:,1),grid_date(:,2),1);
        % log dates of mooring values (with month and year)
        moor_date = datevec(mooring.(moor_nums{m}).time);
        moor_date = datenum(moor_date(:,1),moor_date(:,2),1);
        % first, get index of grid values for mooring dates
        idx_grid = grid_date >= min(moor_date) & grid_date <= max(moor_date);
        % then, get index of mooring values for grid dates
        idx_moor = moor_date >= min(grid_date) & moor_date <= max(grid_date);
        % extract mooring residuals
        mooring_resids = mooring.(moor_nums{m}).pCO2_resids(idx_moor);
        if d == 1
            % for LME-RFR grid
            % calculate differences
            diffs = mooring.(moor_nums{m}).pCO2(idx_moor,1)-...
                mooring.(moor_nums{m}).US_LME_RFR_pCO2_full(idx_grid);
            % date time difference
            dt = mooring.(moor_nums{m}).US_LME_RFR_datetime - datenum(1998,1,15);
            % residuals
            mooring.(moor_nums{m}).US_LME_RFR_pCO2_full_resids = ...
                mooring.(moor_nums{m}).US_LME_RFR_pCO2_full - ...
                (x(1) + x(2).*dt + ...
                x(3).*cos((2.*pi/(365.25/2)).*dt) + x(4).*sin((2.*pi/(365.25/2)).*dt) + ...
                x(5).*cos((2.*pi/(365.25)).*dt) + x(6).*sin((2.*pi/(365.25)).*dt));
            grid_resids = mooring.(moor_nums{m}).US_LME_RFR_pCO2_full_resids(idx_grid);
            % calculate amplitude
            if any(~isnan(mooring.(moor_nums{m}).US_LME_RFR_pCO2_full))
                [~,~,x] = leastsq2(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
                    mooring.(moor_nums{m}).US_LME_RFR_pCO2_full,datenum(1998,1,15),...
                    2,[365.25/2 365.25]);
                amp = sqrt(x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2);
            else
                amp = NaN;
            end
        elseif d == 2
            % for LME-RFR grid (no moorings)
            % calculate differences
            diffs = mooring.(moor_nums{m}).pCO2(idx_moor,1)-...
                mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings(idx_grid);
            % date time difference
            dt = mooring.(moor_nums{m}).US_LME_RFR_datetime - datenum(1998,1,15);
            % residuals
            mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings_resids = ...
                mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings - ...
                (x(1) + x(2).*dt + ...
                x(3).*cos((2.*pi/(365.25/2)).*dt) + x(4).*sin((2.*pi/(365.25/2)).*dt) + ...
                x(5).*cos((2.*pi/(365.25)).*dt) + x(6).*sin((2.*pi/(365.25)).*dt));
                grid_resids = mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings_resids(idx_grid);
            % calculate amplitude
            if any(~isnan(mooring.(moor_nums{m}).US_LME_RFR_pCO2_full))
                [~,~,x] = leastsq2(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
                    mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings,datenum(1998,1,15),...
                    2,[365.25/2 365.25]);
                amp = sqrt(x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2);
            else
                amp = NaN;
            end
        else
            % for other gridded products
            % calculate differences
            diffs = mooring.(moor_nums{m}).pCO2(idx_moor,1)-...
                mooring.(moor_nums{m}).([datasets{d} '_pCO2'])(idx_grid);
            % date time difference
            dt = mooring.(moor_nums{m}).([datasets{d} '_datetime']) - datenum(1998,1,15);
            % residuals
            mooring.(moor_nums{m}).([datasets{d} '_pCO2_resids']) = ...
                mooring.(moor_nums{m}).([datasets{d} '_pCO2']) - ...
                (x(1) + x(2).*dt + ...
                x(3).*cos((2.*pi/(365.25/2)).*dt) + x(4).*sin((2.*pi/(365.25/2)).*dt) + ...
                x(5).*cos((2.*pi/(365.25)).*dt) + x(6).*sin((2.*pi/(365.25)).*dt));
            grid_resids = mooring.(moor_nums{m}).([datasets{d} '_pCO2_resids'])(idx_grid);
            % calculate amplitude
            if any(~isnan(mooring.(moor_nums{m}).US_LME_RFR_pCO2_full))
                [~,~,x] = leastsq2(mooring.(moor_nums{m}).([datasets{d} '_datetime']),...
                    mooring.(moor_nums{m}).([datasets{d} '_pCO2']),datenum(1998,1,15),...
                    2,[365.25/2 365.25]);
                amp = sqrt(x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2);
            else
                amp = NaN;
            end
        end
        % log residual correlation
        idx = ~isnan(mooring_resids) & ~isnan(grid_resids);
        if sum(idx) > 0
            resid_corr(d,m) = corr(mooring_resids(idx),grid_resids(idx));
        else
            resid_corr(d,m) = NaN;
        end
        % log other elements
        mean_err(d,m) = mean(diffs,'omitnan');
        std_err(d,m) = std(diffs,[],'omitnan');
        med_err(d,m) = median(diffs,'omitnan');
        rmse_err(d,m) = sqrt(mean(diffs.^2,'omitnan'));
        amp_diff(d,m) = moor_amp - amp;
        if any(~isnan(diffs))
        iqr_err(d,m) = iqr(diffs);
        else
        iqr_err(d,m) = NaN;
        end
    end
end
mean_table = array2table(mean_err,'VariableNames',moor_nums,'RowNames',datasets);
std_table = array2table(std_err,'VariableNames',moor_nums,'RowNames',datasets);
med_table = array2table(med_err,'VariableNames',moor_nums,'RowNames',datasets);
rmse_table = array2table(rmse_err,'VariableNames',moor_nums,'RowNames',datasets);
iqr_table = array2table(iqr_err,'VariableNames',moor_nums,'RowNames',datasets);
resid_corr_table = array2table(resid_corr,'VariableNames',moor_nums,'RowNames',datasets);
amp_diff_table = array2table(amp_diff,'VariableNames',moor_nums,'RowNames',datasets);
writetable(mean_table,'IndsAndStats/mooring_mean_diffs_all_seaflux.csv');
writetable(std_table,'IndsAndStats/mooring_std_diffs_all_seaflux.csv');
writetable(med_table,'IndsAndStats/mooring_med_diffs_all_seaflux.csv');
writetable(rmse_table,'IndsAndStats/mooring_rmse_diffs_all_seaflux.csv');
writetable(iqr_table,'IndsAndStats/mooring_iqr_diffs_all_seaflux.csv');
writetable(resid_corr_table,'IndsAndStats/mooring_resid_corr_all_seaflux.csv');
writetable(amp_diff_table,'IndsAndStats/mooring_amp_diff_all_seaflux.csv');

%% plot mooring timeseries residuals comparisons
clrs = cbrewer('qual','Set1',8);
datasets = {'US_LME_RFR' 'US_LME_RFR_NM' 'CMEMS_LSCE' 'NIES_FNN' ...
    'MPI_SOMFFN' 'JENA_MLS' 'CSIR_ML6' 'JMA_MLR'};
for m = 1:length(moor_nums)
    % create plot
    f=figure; hold on;
    f.Position(3) = 2*f.Position(3);
    title([num2str(mooring.(moor_nums{m}).US_LME_RFR_lat) char(176) 'N, ' ...
        num2str(mooring.(moor_nums{m}).US_LME_RFR_lon) char(176) 'E']);
    p.p1=plot(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).US_LME_RFR_pCO2_full_resids,...
        '-','color',clrs(1,:),'linewidth',2);
    p.p2=plot(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).US_LME_RFR_pCO2_no_moorings_resids,...
        '--','color',clrs(2,:),'linewidth',2);
    for d = 3:8
        p.(['p' num2str(d)])=...
            plot(mooring.(moor_nums{m}).([datasets{d} '_datetime']),...
            mooring.(moor_nums{m}).([datasets{d} '_pCO2_resids']),...
            '-','color',[clrs(d,:) 0.4],'linewidth',1);
    end
    s1=scatter(mooring.(moor_nums{m}).US_LME_RFR_datetime,...
        mooring.(moor_nums{m}).pCO2_resids,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    % add figure properties
    xlim([min(mooring.(moor_nums{m}).US_LME_RFR_datetime) ...
        max(mooring.(moor_nums{m}).US_LME_RFR_datetime)]);
    ylabel('{\itp}CO_{2} (\muatm)');
    datetick('x','yyyy','keeplimits');
    if mooring.(moor_nums{m}).lat == 59.875
        ylim([150 500]);
    end
    legend([p.p1 p.p2 p.p3 p.p4 p.p5 p.p6 p.p7 p.p8],...
        {'LME-RFR' 'LME-RFR_{NM}' 'CMEMS-LSCE' 'NIES-FNN' ...
        'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
        'Location','northwest','NumColumns',4);
    % export figure
    % exportgraphics(gcf,['Figures/moorings/m' moor_nums{m} '_resid.png']);
    export_fig(gcf,['Figures/moorings/m' moor_nums{m} '_resid.png'],'-transparent');
    close
end

%% plot box and whisker figures
% index to only within LMEs
idx = ~isnan(mean_err(1,:));
% establish figure
figure;
set(gcf,'Position',[100 100 1000 800])
gap = [0.09,0.08];
% median errors
subtightplot(2,2,1,gap);
h=boxplot(med_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-100 100]);
ylabel('Median \Delta{\itp}CO_{2} (\muatm)');
% IQR errors
subtightplot(2,2,2,gap);
h=boxplot(iqr_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
ylim([0 140]);
ylabel('\Delta{\itp}CO_{2} IQR (\muatm)');
% amplitude difference errors
subtightplot(2,2,3,gap);
h=boxplot(amp_diff(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-20 130]);
ylabel('\Delta {\itp}CO_{2} Amplitude (\muatm)');
% residual correlation errors
subtightplot(2,2,4,gap)
h=boxplot(resid_corr(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-0.5 1]);
ylabel('{\itp}CO_{2} Residual Correlation');
% export plots
% exportgraphics(gcf,'Figures/boxplots_moorings.png');
export_fig(gcf,'Figures/boxplots_moorings.png','-transparent');
close;

%% plot box and whisker figure
% index to only within LMEs
idx = ~isnan(mean_err(1,:));
% median errors
h=boxplot(med_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-100 100]);
ylabel('Median \Delta{\itp}CO_{2} (\muatm)');
% exportgraphics(gcf,['Figures/median_boxplot_moorings.png']);
export_fig(gcf,['Figures/median_boxplot_moorings.png'],'-transparent');
close;
% mean errors
h=boxplot(mean_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-100 100]);
ylabel('Mean \Delta{\itp}CO_{2} (\muatm)');
% exportgraphics(gcf,['Figures/mean_boxplot_moorings.png']);
export_fig(gcf,['Figures/mean_boxplot_moorings.png'],'-transparent');
close;
% RMSE errors
figure;
h=boxplot(rmse_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
ylim([0 140]);
ylabel('\Delta{\itp}CO_{2} RMSD (\muatm)');
% exportgraphics(gcf,['Figures/rmse_boxplot_moorings.png']);
export_fig(gcf,['Figures/rmse_boxplot_moorings.png'],'-transparent');
close;
% stdev errors
figure;
h=boxplot(std_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
ylim([0 140]);
ylabel('{\itp}CO_{2} StDev (\muatm)');
% exportgraphics(gcf,['Figures/std_boxplot_moorings.png']);
export_fig(gcf,['Figures/std_boxplot_moorings.png'],'-transparent');
close;
% IQR errors
figure; set(gcf,'Position',[100 100 700 600]);
h=boxplot(iqr_err(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
ylim([0 140]);
ylabel('\Delta{\itp}CO_{2} IQR (\muatm)');
% exportgraphics(gcf,['Figures/iqr_boxplot_moorings.png']);
export_fig(gcf,['Figures/iqr_boxplot_moorings.png'],'-transparent');
close;
% residual correlation errors
figure;
h=boxplot(resid_corr(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-0.5 1]);
ylabel('{\itp}CO_{2} Residual Correlation');
% exportgraphics(gcf,['Figures/resid_corr_boxplot_moorings.png']);
export_fig(gcf,['Figures/resid_corr_boxplot_moorings.png'],'-transparent');
close;
% amplitude difference errors
figure;
h=boxplot(amp_diff(:,idx)',{'RFR-LME' 'RFR-LME-NM' 'CMEMS-LSCE' ...
    'NIES-FNN' 'MPI-SOMFFN' 'JENA-MLS' 'CSIR-ML6' 'JMA-MLR'},...
    'BoxStyle','outline','Symbol','.k','Colors','k');
set(h,'linew',1.5);
yline(0);
ylim([-20 130]);
ylabel('\Delta {\itp}CO_{2} Amplitude (\muatm)');
% exportgraphics(gcf,['Figures/amp_diff_boxplot_moorings.png']);
export_fig(gcf,['Figures/amp_diff_boxplot_moorings.png'],'-transparent');
close;

%% plot mean map with moorings
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(flipud(slanCM('romao')));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)';
cbarrow;
% plot regions
for n = 1:length(region)
    vars_grid = load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');
    z = mean(vars_grid.OAI_grid.(region{n}).pCO2,3,'omitnan')';
    contourfm(vars_grid.OAI_grid.(region{n}).lat,vars_grid.OAI_grid.(region{n}).lon,...
        z,295:(475-295)/200:475,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
plot_lme_borders(region,lme_shape,lme_idx);
% plot land
plot_land('map');
mlabel off
% plot pCO2 from moorings
for m = 1:length(moor_nums)
    scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
        40,mean(mooring.(moor_nums{m}).pCO2,'omitnan'),'filled',...
        'MarkerEdgeColor','k')
end
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
% exportgraphics(gcf,['Figures/full/pCO2_with_moorings.png']);
export_fig(gcf,['Figures/full/pCO2_with_moorings.png'],'-transparent');
close

%% plot trend map with moorings
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
caxis([-3 3]);
colormap(cmocean('balance','pivot',0));
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itp}CO_{2(RFR-LME)} Trend (\muatm yr^{-1})';
cbarrow;
% plot regions
for n = 1:length(region)
    vars_grid = load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');
    tr = nan(vars_grid.OAI_grid.(region{n}).dim.x,vars_grid.OAI_grid.(region{n}).dim.y);
    for a = 1:vars_grid.OAI_grid.(region{n}).dim.x
        for b = 1:vars_grid.OAI_grid.(region{n}).dim.y
            if sum(~isnan(squeeze(vars_grid.OAI_grid.(region{n}).pCO2(a,b,:))))>200 % if more than 200 months with observations
            [~,~,x,~] = leastsq2(vars_grid.OAI_grid.(region{n}).month,...
                squeeze(vars_grid.OAI_grid.(region{n}).pCO2(a,b,:)),0,2,[6 12]);
            tr(a,b) = x(2)*12;
            else
                tr(a,b) = NaN;
            end
        end
    end
    contourfm(vars_grid.OAI_grid.(region{n}).lat,vars_grid.OAI_grid.(region{n}).lon,...
        tr',-3:6/200:3,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
plot_lme_borders(region,lme_shape,lme_idx);
% plot land
plot_land('map');
mlabel off
% plot pCO2 from moorings
for m = 1:length(moor_nums)
    [~,~,x] = leastsq2(mooring.(moor_nums{m}).time,...
        mooring.(moor_nums{m}).pCO2,mooring.(moor_nums{m}).time(1),2,[365/2 365]);
    mooring.(moor_nums{m}).trend = x(2)*365;
    scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
        40,mooring.(moor_nums{m}).trend,'filled','MarkerEdgeColor','k');
end
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
% exportgraphics(gcf,['Figures/full/pCO2_tr_with_moorings.png']);
export_fig(gcf,['Figures/full/pCO2_tr_with_moorings.png'],'-transparent');
close

%% caluclate climatologies
% for mnth = 1:12
%     mooring_clim.(['m' num2str(m)]).pCO2(mnth) = ...
%         mean(mooring.(['m' num2str(m)]).pCO2(mnth:12:end),'omitnan');
%     no_moorings_clim(mnth) = ...
%         mean(US_LME_RFR.pCO2_no_moorings(a(m),b(m),mnth:12:end),3,'omitnan');
%     moorings_clim(mnth) = ...
%         mean(US_LME_RFR.pCO2_full(a(m),b(m),mnth:12:end),3,'omitnan');
% end
% 
% % add stats
% del_no_moorings = mooring.(['m' num2str(m)]).pCO2 - ...
%     squeeze(US_LME_RFR.pCO2_no_moorings(a(m),b(m),:));
% del_full = mooring.(['m' num2str(m)]).pCO2 - ...
%     squeeze(US_LME_RFR.pCO2_full(a(m),b(m),:));
% avg_del_no_moorings = mean(del_no_moorings,'omitnan');
% avg_del_full = mean(del_full,'omitnan');
% std_del_no_moorings = std(del_no_moorings,[],'omitnan');
% std_del_full = std(del_full,[],'omitnan');
% legend([p1 p2],{['No Moorings: ' num2str(round(avg_del_no_moorings,1)) ' ' char(177) ' ' ...
%     num2str(round(std_del_no_moorings,1))],['Moorings: ' ...
%     num2str(round(avg_del_full,1)) ' ' char(177) ' ' num2str(round(std_del_full,1))]},...
%     'location','northwest','fontsize',12);
% 
% % display mooring number
% disp(['Moor. # ' num2str(m)]);
% % display means
% disp(['Moor. Mean = ' num2str(mean(mooring_clim.(['m' num2str(m) '_pCO2']),'omitnan'))]);
% disp(['RFR Mean (no moor.) = ' num2str(mean(no_moorings_clim,'omitnan'))]);
% disp(['RFR Mean (w/ moor.) = ' num2str(mean(moorings_clim,'omitnan'))]);
% % display amplitudes
% disp(['Moor. Amp. = ' num2str(max(mooring_clim.(['m' num2str(m) '_pCO2']))-...
%     min(mooring_clim.(['m' num2str(m) '_pCO2'])))]);
% disp(['RFR Amp. (no moor.) = ' num2str(max(no_moorings_clim)-...
%     min(no_moorings_clim))]);
% disp(['RFR Amp. (w/ moor.) = ' num2str(max(moorings_clim)-...
%     min(moorings_clim))]);