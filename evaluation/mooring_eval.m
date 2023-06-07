% Evaluate US LME OA indicators against mooring observations

%% Import mooring data from csv files
ImportMoorings

%% load US LME data
date = '20-May-2023';
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

%% extract time series from RFR data
for n = 1:length(moornames)
    MOORING.(moornames{n}).datetime = datenum(MOORING.(moornames{n}).date);
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    lon_idx = find(abs(LME_RFR.Lon - moor_lon) == min(abs(LME_RFR.Lon - moor_lon)));
    lat_idx = find(abs(LME_RFR.Lat - moor_lat) == min(abs(LME_RFR.Lat - moor_lat)));
    MOORING.(moornames{n}).LME_RFR_fCO2 = squeeze(LME_RFR.fCO2(lon_idx,lat_idx,:));
    MOORING.(moornames{n}).LME_RFR_ufCO2 = squeeze(LME_RFR_u.ufCO2(lon_idx,lat_idx,:));
    MOORING.(moornames{n}).LME_RFR_datetime = LME_RFR.Time;
end

%% extract time series from CMEMS-LSCE data
for n = 1:length(moornames)
    MOORING.(moornames{n}).datetime = datenum(MOORING.(moornames{n}).date);
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    londiffs = sort(abs(CMEMS_LSCE.longitude - moor_lon));
    latdiffs = sort(abs(CMEMS_LSCE.latitude - moor_lat));
    lon_idx = find(abs(CMEMS_LSCE.longitude - moor_lon) == londiffs(1));
    lat_idx = find(abs(CMEMS_LSCE.latitude - moor_lat) == latdiffs(1));
    % try moving lat a bit
    if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1))
        lat_idx = find(abs(CMEMS_LSCE.latitude - moor_lat) == latdiffs(2));
    end
    % try moving lon a bit
    if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1))
        lon_idx = find(abs(CMEMS_LSCE.longitude - moor_lon) == londiffs(2));
    end
    MOORING.(moornames{n}).CMEMS_LSCE_pCO2 = squeeze(CMEMS_LSCE.spco2(lon_idx,lat_idx,:));
    MOORING.(moornames{n}).CMEMS_LSCE_datetime = CMEMS_LSCE.time;
end

%% extract time series from OceanSODA-ETHZ data
for n = 1:length(moornames)
    MOORING.(moornames{n}).datetime = datenum(MOORING.(moornames{n}).date);
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    londiffs = sort(abs(SODA_ETHZ.lon - moor_lon));
    latdiffs = sort(abs(SODA_ETHZ.lat - moor_lat));
    lon_idx = find(abs(SODA_ETHZ.lon - moor_lon) == londiffs(1));
    lat_idx = find(abs(SODA_ETHZ.lat - moor_lat) == latdiffs(1));
    % try moving lat a bit
    if isnan(SODA_ETHZ.spco2(lon_idx,lat_idx,1))
        lat_idx = find(abs(SODA_ETHZ.lat - moor_lat) == latdiffs(2));
    end
    % try moving lon a bit
    if isnan(SODA_ETHZ.spco2(lon_idx,lat_idx,1))
        lon_idx = find(abs(SODA_ETHZ.lon - moor_lon) == londiffs(2));
    end
    MOORING.(moornames{n}).SODA_ETHZ_pCO2 = squeeze(SODA_ETHZ.spco2(lon_idx,lat_idx,:));
    MOORING.(moornames{n}).SODA_ETHZ_datetime = SODA_ETHZ.time;
end

%% extract time series from MPI-SOMFFN data
for n = 1:length(moornames)
    MOORING.(moornames{n}).datetime = datenum(MOORING.(moornames{n}).date);
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    londiffs = sort(abs(MPI_SOMFFN.lon - moor_lon));
    latdiffs = sort(abs(MPI_SOMFFN.lat - moor_lat));
    lon_idx = find(abs(MPI_SOMFFN.lon - moor_lon) == londiffs(1));
    lat_idx = find(abs(MPI_SOMFFN.lat - moor_lat) == latdiffs(1));
    % try moving lat a bit
    if isnan(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,1))
        lat_idx = find(abs(MPI_SOMFFN.lat - moor_lat) == latdiffs(2));
    end
    % try moving lon a bit
    if isnan(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,1))
        lon_idx = find(abs(MPI_SOMFFN.lon - moor_lon) == londiffs(2));
    end
    MOORING.(moornames{n}).MPI_SOMFFN_pCO2 = squeeze(MPI_SOMFFN.spco2_raw(lon_idx,lat_idx,:));
    MOORING.(moornames{n}).MPI_SOMFFN_datetime = MPI_SOMFFN.time;
end

%% plot mooring against RFR
for n = 1:length(moornames)
    % plot fCO2 time series
    clrs = cbrewer('qual','Set1',4);
    f=figure; hold on;
    f.Position(3) = 2*f.Position(3);
    p1=plot(MOORING.(moornames{n}).LME_RFR_datetime,MOORING.(moornames{n}).LME_RFR_fCO2,...
        'color',clrs(1,:),'linewidth',2);
    fill([MOORING.(moornames{n}).LME_RFR_datetime;flipud(MOORING.(moornames{n}).LME_RFR_datetime)],...
        [MOORING.(moornames{n}).LME_RFR_fCO2+MOORING.(moornames{n}).LME_RFR_ufCO2;...
        flipud(MOORING.(moornames{n}).LME_RFR_fCO2-MOORING.(moornames{n}).LME_RFR_ufCO2)],...
        clrs(1,:),'facealpha',0.25,'linestyle','none')
%     p2=plot(MOORING.(moornames{n}).CMEMS_LSCE_datetime,MOORING.(moornames{n}).CMEMS_LSCE_pCO2,...
%         'color',clrs(2,:),'linewidth',2);
%     p3=plot(MOORING.(moornames{n}).SODA_ETHZ_datetime,MOORING.(moornames{n}).SODA_ETHZ_pCO2,...
%         'color',clrs(3,:),'linewidth',2);
%     p4=plot(MOORING.(moornames{n}).MPI_SOMFFN_datetime,MOORING.(moornames{n}).MPI_SOMFFN_pCO2,...
%         'color',clrs(4,:),'linewidth',2);
    s1=scatter(MOORING.(moornames{n}).datetime_monthly,MOORING.(moornames{n}).fCO2SW_monthly(:,1),'ko','filled');
%     xlim([min(MOORING.(moornames{n}).LME_RFR_datetime) max(MOORING.(moornames{n}).LME_RFR_datetime)]);
    xlim([MOORING.(moornames{n}).LME_RFR_datetime(96) max(MOORING.(moornames{n}).LME_RFR_datetime)]);
    datetick('x','yyyy','keeplimits');
%     legend([p1 p2 p3 p4 s1],{'LME-RFR' 'CMEMS-LSCE' 'OceanSODA-ETHZ' 'MPI-SOMFFN' 'Mooring'},...
%         'Location','north','NumColumns',5);
    legend([p1 s1],{'LME-RFR' 'Mooring'},...
        'Location','north','NumColumns',2);
    exportgraphics(f,['Figures/mooring_comp_' moornames{n} '.png']);
    close
    % plot absolute difference between mooring and RFR
end

%% plot mooring averages on map
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([28 62],[194 244]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(parula(18));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itf}CO_{2} (\muatm)';
cbarrow;
% plot land
bordersm('alaska','facecolor',rgb('gray'))
bordersm('continental us','facecolor',rgb('gray'))
bordersm('canada','facecolor',rgb('light grey'))
mlabel off
% plot regions
for n = 1:4
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    z = mean(OAI_grid.(region{n}).fCO2,3,'omitnan')';
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        z,295:10:475,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
for n = 1:4
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add Moorings to plot
for n = 1:length(moornames)
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    moor_fCO2 = mean(MOORING.(moornames{n}).fCO2SW_clim(:,1),'omitnan');
    scatterm(moor_lat,moor_lon,100,moor_fCO2,'filled','markeredgecolor','k');
end
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Mooring_eval_map.png');
% print means
for n = 1:length(moornames)
    % Mooring mean
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    moor_time_min = min(MOORING.(moornames{n}).datetime);
    moor_time_max = max(MOORING.(moornames{n}).datetime);
    disp(['Mooring Mean (' moornames{n} ') = ' ...
        num2str(round(MOORING.(moornames{n}).fCO2SW_annmean,1)) ' uatm']);
    % RFR-LME mean
    lon_idx = find(abs(LME_RFR.Lon - moor_lon) == min(abs(LME_RFR.Lon - moor_lon)));
    lat_idx = find(abs(LME_RFR.Lat - moor_lat) == min(abs(LME_RFR.Lat - moor_lat)));
    t1_idx = find(abs(LME_RFR.Time - moor_time_min) == min(abs(LME_RFR.Time - moor_time_min)));
    t2_idx = find(abs(LME_RFR.Time - moor_time_max) == min(abs(LME_RFR.Time - moor_time_max)));
    disp(['RFR-LME Mean (' moornames{n} ') = ' ...
        num2str(round(mean(LME_RFR.fCO2(lon_idx,lat_idx,t1_idx:t2_idx),3,'omitnan'),1)) ' uatm']);
    % CMEMS-LSCE mean
    londiffs = sort(abs(CMEMS_LSCE.longitude - moor_lon));
    latdiffs = sort(abs(CMEMS_LSCE.latitude - moor_lat));
    lon_idx = find(abs(CMEMS_LSCE.longitude - moor_lon) == londiffs(1));
    lat_idx = find(abs(CMEMS_LSCE.latitude - moor_lat) == latdiffs(1));
    t1_idx = find(abs(CMEMS_LSCE.time - moor_time_min) == min(abs(CMEMS_LSCE.time - moor_time_min)));
    t2_idx = find(abs(CMEMS_LSCE.time - moor_time_max) == min(abs(CMEMS_LSCE.time - moor_time_max)));
    if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1)) % try moving lat a bit
        lat_idx = find(abs(CMEMS_LSCE.latitude - moor_lat) == latdiffs(2));
    end
    if isnan(CMEMS_LSCE.spco2(lon_idx,lat_idx,1)) % try moving lon a bit
        lon_idx = find(abs(CMEMS_LSCE.longitude - moor_lon) == londiffs(2));
    end
    disp(['CMEMS-LSCE Mean (' moornames{n} ') = ' ...
        num2str(round(mean(CMEMS_LSCE.spco2(lon_idx,lat_idx,t1_idx:t2_idx),3,'omitnan'),1)) ' uatm']);
end

%% plot mooring amplitudes on map
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([28 62],[194 244]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(cmocean('thermal',15));
caxis([1 2.5]);
c.TickLength = 0;
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.String = 'Sea Surface {\itf}CO_{2} Amplitude (\muatm)';
% plot land
bordersm('alaska','facecolor',rgb('gray'))
bordersm('continental us','facecolor',rgb('gray'))
bordersm('canada','facecolor',rgb('light grey'))
mlabel off
% plot regions
for n = 1:4
    % load
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % calculate climatology
    OAI_grid.(region{n}).fCO2_clim = ...
        nan(OAI_grid.(region{n}).dim.x,OAI_grid.(region{n}).dim.y,12);
    for m = 1:12
        OAI_grid.(region{n}).fCO2_clim(:,:,m) = ...
            mean(OAI_grid.(region{n}).fCO2(:,:,m:12:end),3,'omitnan');
    end
    % calculate amplitude
    OAI_grid.(region{n}).fCO2_amp = ...
        max(OAI_grid.(region{n}).fCO2_clim,[],3) - ...
        min(OAI_grid.(region{n}).fCO2_clim,[],3);
    % plot
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        log10(OAI_grid.(region{n}).fCO2_amp)',1:0.1:2.5,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
for n = 1:4
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add Moorings to plot
for n = 1:length(moornames)
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    scatterm(moor_lat,moor_lon,100,log10(MOORING.(moornames{n}).fCO2SW_amplitude),...
        'filled','markeredgecolor','k');
end
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Mooring_eval_amp_map.png');
% print amplitudes
for n = 1:length(moornames)
    % Mooring amp
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    moor_time_min = min(MOORING.(moornames{n}).datetime);
    moor_time_max = max(MOORING.(moornames{n}).datetime);
    disp(['Mooring Amp. (' moornames{n} ') = ' ...
        num2str(round(MOORING.(moornames{n}).fCO2SW_amplitude,1)) ' uatm']);
    % RFR-LME climatology
    LME_RFR.fCO2_clim = ...
        nan(length(LME_RFR.Lon),length(LME_RFR.Lat),12);
    for m = 1:12
        LME_RFR.fCO2_clim(:,:,m) = ...
            mean(LME_RFR.fCO2(:,:,m:12:end),3,'omitnan');
    end
    % RFR-LME amplitude
    LME_RFR.fCO2_amp = ...
        max(LME_RFR.fCO2_clim,[],3) - ...
        min(LME_RFR.fCO2_clim,[],3);
    lon_idx = find(abs(LME_RFR.Lon - moor_lon) == min(abs(LME_RFR.Lon - moor_lon)));
    lat_idx = find(abs(LME_RFR.Lat - moor_lat) == min(abs(LME_RFR.Lat - moor_lat)));
    t1_idx = find(abs(LME_RFR.Time - moor_time_min) == min(abs(LME_RFR.Time - moor_time_min)));
    t2_idx = find(abs(LME_RFR.Time - moor_time_max) == min(abs(LME_RFR.Time - moor_time_max)));
    disp(['RFR-LME Amp. (' moornames{n} ') = ' ...
        num2str(round(mean(LME_RFR.fCO2_amp(lon_idx,lat_idx),3,'omitnan'),1)) ' uatm']);
end

%% plot mooring trends on map
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([28 62],[194 244]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
clim([-1 4]);
colormap(cmocean('balance',20,'pivot',0));
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itf}CO_{2} Trend (\muatm yr^{-1})';
% plot land
bordersm('alaska','facecolor',rgb('gray'))
bordersm('continental us','facecolor',rgb('gray'))
bordersm('canada','facecolor',rgb('light grey'))
mlabel off
% plot regions
for n = 1:4
    % load
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % calculate trends
    OAI_grid.(region{n}).fCO2_trend = ...
        nan(OAI_grid.(region{n}).dim.x,OAI_grid.(region{n}).dim.y);
    for a = 1:length(OAI_grid.(region{n}).lon)
        for b = 1:length(OAI_grid.(region{n}).lat)
            mod = polyfit(OAI_grid.(region{n}).month,...
                squeeze(OAI_grid.(region{n}).fCO2(a,b,:)),1);
            OAI_grid.(region{n}).fCO2_trend(a,b) = mod(1)*12;
        end
    end
    % plot
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).fCO2_trend',-1:0.25:4,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
for n = 1:4
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add Moorings to plot
for n = 1:length(moornames)
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    scatterm(moor_lat,moor_lon,100,MOORING.(moornames{n}).fCO2SW_trend,...
        'filled','markeredgecolor','k');
end
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Mooring_eval_trend_map.png');
% print amplitudes
for n = 1:length(moornames)
    % Mooring trend
    moor_lon = mean(MOORING.(moornames{n}).lon,'omitnan');
    moor_lat = mean(MOORING.(moornames{n}).lat,'omitnan');
    disp(['Mooring Trend (' moornames{n} ') = ' ...
        num2str(round(MOORING.(moornames{n}).fCO2SW_trend,1)) ' uatm/yr']);
    % RFR-LME trend
    lon_idx = find(abs(LME_RFR.Lon - moor_lon) == min(abs(LME_RFR.Lon - moor_lon)));
    lat_idx = find(abs(LME_RFR.Lat - moor_lat) == min(abs(LME_RFR.Lat - moor_lat)));
    mod = polyfit(LME_RFR.Time-LME_RFR.Time(1),squeeze(LME_RFR.fCO2(lon_idx,lat_idx,:)),1);
    disp(['RFR-LME Trend (' moornames{n} ') = ' ...
        num2str(round(mean(mod(1)*365,3,'omitnan'),1)) ' uatm/yr']);
end
