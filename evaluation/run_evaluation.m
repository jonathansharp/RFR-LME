% Run scripts to evaluate US LME OA indicators

addpath(genpath(pwd));

% RFR-LME file data!
date = '26-Feb-2024';
%date = '06-Oct-2023';

new = 1;

%% load CODAP data
load('CODAP_NA/CODAP_NA_v2020_G2format.mat')
idx = C1.pressure <= 10 & C1.talkf == 2 & C1.tco2f == 2 & ...
    C1.salinityf == 2 & C1.silicatef == 2 & C1.phosphatef == 2;
CODAP.lat = C1.latitude(idx);
CODAP.lon = C1.longitude(idx);
CODAP.lon = convert_lon(CODAP.lon);
CODAP.time = datenum(C1.year(idx),C1.month(idx),C1.day(idx));
CODAP.prs = C1.pressure(idx);
CODAP.sal = C1.salinity(idx);
CODAP.tmp = C1.temperature(idx);
CODAP.sil = C1.silicate(idx);
CODAP.phos = C1.phosphate(idx);
CODAP.DIC = C1.tco2(idx);
CODAP.TA = C1.talk(idx);
clear C1 idx Headers Units

%% load GLODAP data
load('GLODAPv2.2022/GLODAPv2.2022_Merged_Master_File.mat')
idx = G2pressure <= 10 & G2talkf == 2 & G2tco2f == 2 & ...
    G2salinityf == 2 & G2silicatef == 2 & G2phosphatef == 2;
GLODAP.lat = G2latitude(idx);
GLODAP.lon = G2longitude(idx);
GLODAP.lon = convert_lon(GLODAP.lon);
GLODAP.time = datenum(G2year(idx),G2month(idx),G2day(idx));
GLODAP.prs = G2pressure(idx);
GLODAP.sal = G2salinity(idx);
GLODAP.tmp = G2temperature(idx);
GLODAP.sil = G2silicate(idx);
GLODAP.phos = G2phosphate(idx);
GLODAP.DIC = G2tco2(idx);
GLODAP.TA = G2talk(idx);
clear idx G2* expocode expocodeno

%% remove CODAP data outside LMEs
if new == 1
    RFR_LME.time = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'time');
else
    RFR_LME.time = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Time');
end
% determine LME index
define_regions_eiwg
idx_tmp = nan(length(CODAP.lat),length(region));
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_tmp_CODAP(:,n) = inpolygon(CODAP.lon,CODAP.lat,tmp_lon,tmp_lat);
    idx_tmp_GLODAP(:,n) = inpolygon(GLODAP.lon,GLODAP.lat,tmp_lon,tmp_lat);
end
idx_CODAP = any(idx_tmp_CODAP,2);
idx_GLODAP = any(idx_tmp_GLODAP,2);
% also remove data outside time limits
idx_CODAP = idx_CODAP & CODAP.time > min(RFR_LME.time) & CODAP.time < max(RFR_LME.time);
idx_GLODAP = idx_GLODAP & GLODAP.time > min(RFR_LME.time) & GLODAP.time < max(RFR_LME.time);
% remove outside data (CODAP)
vars = fieldnames(CODAP);
for v = 1:length(vars)
    CODAP.(vars{v}) = CODAP.(vars{v})(idx_CODAP);
end
% remove outside data (GLODAP)
vars = fieldnames(GLODAP);
for v = 1:length(vars)
    GLODAP.(vars{v}) = GLODAP.(vars{v})(idx_GLODAP);
end
clear idx_CODAP idx_tmp_CODAP idx_GLODAP idx_tmp_GLODAP vars v tmp_lon tmp_lat

%% calculate fCO2, H, and CO3 (CODAP)
carb = CO2SYS(CODAP.TA,CODAP.DIC,1,2,CODAP.sal,CODAP.tmp,NaN,CODAP.prs,...
    NaN,CODAP.sil,CODAP.phos,0,0,1,10,1,2,2);
carb_err = errors(CODAP.TA,CODAP.DIC,1,2,CODAP.sal,CODAP.tmp,NaN,CODAP.prs,...
    NaN,CODAP.sil,CODAP.phos,0,0,2,2,0.01,0.01,0.1,0.01,0,0,'','',0,1,10,1,2,2);
carb(carb==-999)=NaN;
carb_err(carb_err==-999)=NaN;
CODAP.pco2 = carb(:,4);
CODAP.pco2_e = carb_err(:,4);
CODAP.fco2 = carb(:,5);
CODAP.fco2_e = carb_err(:,5);
CODAP.pH = carb(:,3);
CODAP.pH_e = -log10(carb(:,15)./10^6) -  (-log10((carb(:,15)./10^6+carb_err(:,3)./10^9)));
CODAP.H = carb(:,15).*10^3;
CODAP.H_e = carb_err(:,3).*10^3;
CODAP.CO3 = carb(:,7);
CODAP.CO3_e = carb_err(:,7);
CODAP.OmA = carb(:,18);
CODAP.OmA_e = carb_err(:,11);
CODAP.OmC = carb(:,17);
CODAP.OmC_e = carb_err(:,10);
CODAP.RF = carb(:,16);
CODAP.RF_e = carb_err(:,9);
clear carb carb_err

%% calculate fCO2, H, CO3, Omega, and RF (GLODAP)
carb = CO2SYS(GLODAP.TA,GLODAP.DIC,1,2,GLODAP.sal,GLODAP.tmp,NaN,GLODAP.prs,...
    NaN,GLODAP.sil,GLODAP.phos,0,0,1,10,1,2,2);
carb_err = errors(GLODAP.TA,GLODAP.DIC,1,2,GLODAP.sal,GLODAP.tmp,NaN,GLODAP.prs,...
    NaN,GLODAP.sil,GLODAP.phos,0,0,2,2,0,0,0,0,0,0,'','',0,1,10,1,2,2);
carb(carb==-999)=NaN;
carb_err(carb_err==-999)=NaN;
GLODAP.pco2 = carb(:,4);
GLODAP.pco2_e = carb_err(:,4);
GLODAP.fco2 = carb(:,5);
GLODAP.fco2_e = carb_err(:,5);
GLODAP.pH = carb(:,3);
GLODAP.pH_e = -log10(carb(:,15)./10^6) -  (-log10((carb(:,15)./10^6+carb_err(:,3)./10^9)));
GLODAP.H = carb(:,15).*10^3;
GLODAP.H_e = carb_err(:,3).*10^3;
GLODAP.CO3 = carb(:,7);
GLODAP.CO3_e = carb_err(:,7);
GLODAP.OmA = carb(:,18);
GLODAP.OmA_e = carb_err(:,11);
GLODAP.OmC = carb(:,17);
GLODAP.OmC_e = carb_err(:,20);
GLODAP.RF = carb(:,16);
GLODAP.RF_e = carb_err(:,9);
clear carb carb_err

%% load US LME data
if new == 1
% New NetCDF files
RFR_LME.lat = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'lat');
RFR_LME.lon = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'lon');
RFR_LME.time = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'time');
RFR_LME.pco2 = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'pco2');
RFR_LME.fco2 = ncread(['Data/NetCDFs_' date '/US_RFR_LME_fCO2.nc'],'fco2');
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
RFR_LME.lat = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Lat');
RFR_LME.lon = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Lon');
RFR_LME.time = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'Time');
RFR_LME.pco2 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'pCO2');
RFR_LME.fco2 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'fCO2');
RFR_LME.TA = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'TA');
RFR_LME.DIC = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'DIC');
RFR_LME.pH = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'pH');
RFR_LME.OmA = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'OmA');
RFR_LME.OmC = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'OmC');
RFR_LME.H = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'H');
RFR_LME.CO3 = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'CO3');
RFR_LME.RF = ncread(['Data/US_RFR_LME_Inds_' date '.nc'],'RF');
end

%% plot data locations
% define_regions_eiwg
% % initialize figure
% figure('visible','on'); box on; hold on;
% worldmap([-18 82],[140 302]);
% setm(gca,'MapProjection','robinson','MLabelParallel','south');
% set(gcf,'position',[100 100 900 600]);
% set(gca,'fontsize',16);
% % figure properties
% c=colorbar('location','southoutside');
% colormap(flipud(slanCM('romao')));
% caxis([0 5]);
% c.TickLength = 0;
% c.Label.String = 'Sea Surface \Omega_{ar(RFR-LME)} (\muatm)';
% cbarrow;
% % plot regions
% for n = 1:length(region)
%      load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
%     z = mean(OAI_grid.(region{n}).OmA,3,'omitnan')';
%     contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
%         z,0:5/200:5,'LineStyle','none');
%     clear vars_grid z
% end
% % plot borders around regions
% for n = 1:length(region)
%     if n <= 11
%         tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
%     else
%         tmp_lon = lme_shape(lme_idx.(region{n})).X';
%     end
%     tmp_lat = lme_shape(lme_idx.(region{n})).Y';
%     plotm(tmp_lat,tmp_lon,'k','linewidth',1);
% end
% % plot land
% plot_land('map');
% mlabel off
% % add CODAP to plot
% scatterm(CODAP.lat,CODAP.lon,10,'magenta','o');
% % add GLODAP to plot
% scatterm(GLODAP.lat,GLODAP.lon,10,'green','.');
% % save figure
% if ~isfolder('Figures'); mkdir('Figures'); end
% exportgraphics(gcf,'Figures/CODAP_GLODAP_eval_map.png');

%% variable information
edges = {1900:10:2300;300:5:500;300:5:500;2000:10:2400;7.8:0.0125:8.3;...
    1:0.1:5;1:0.1:5;7:0.2:14;100:2.5:250;9:0.1:17};
var_type = {'DIC' 'fco2' 'pco2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}_{CO2}' '{\itp}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 1 3 2 2 1 1 2];

%% co-locate each CODAP point with an LME-RFR grid cell
% pre-allocate
for v = 1:length(var_type)
    CODAP.([var_type{v} '_grid']) = nan(size(CODAP.(var_type{v})));
end
% identify and fill gridded values
for n = 1:length(CODAP.lat)
    % lon
    idx_lon = find(abs(RFR_LME.lon - CODAP.lon(n)) == ...
        min(abs(RFR_LME.lon - CODAP.lon(n))));
    % lat
    idx_lat = find(abs(RFR_LME.lat - CODAP.lat(n)) == ...
        min(abs(RFR_LME.lat - CODAP.lat(n))));
    % time
    idx_time = find(abs(RFR_LME.time - CODAP.time(n)) == ...
        min(abs(RFR_LME.time - CODAP.time(n))));
    for v = 1:length(var_type)
        CODAP.([var_type{v} '_grid'])(n) = ...
            RFR_LME.(var_type{v})(idx_lon(1),idx_lat(1),idx_time(1));
    end
end
% calculate differences
for v = 2:length(var_type)
    CODAP.([var_type{v} '_del']) = ...
        CODAP.(var_type{v}) - CODAP.([var_type{v} '_grid']);
end

%% co-locate each GLODAP point with an LME-RFR grid cell
% pre-allocate
for v = 1:length(var_type)
    GLODAP.([var_type{v} '_grid']) = nan(size(GLODAP.(var_type{v})));
end
% identify and fill gridded values
for n = 1:length(GLODAP.lat)
    % lon
    idx_lon = find(abs(RFR_LME.lon - GLODAP.lon(n)) == ...
        min(abs(RFR_LME.lon - GLODAP.lon(n))));
    % lat
    idx_lat = find(abs(RFR_LME.lat - GLODAP.lat(n)) == ...
        min(abs(RFR_LME.lat - GLODAP.lat(n))));
    % time
    idx_time = find(abs(RFR_LME.time - GLODAP.time(n)) == ...
        min(abs(RFR_LME.time - GLODAP.time(n))));
    for v = 1:length(var_type)
        GLODAP.([var_type{v} '_grid'])(n) = ...
            RFR_LME.(var_type{v})(idx_lon(1),idx_lat(1),idx_time(1));
    end
end
% eliminate NaN and -999
for v = 1:length(var_type)
    idx = ~isnan(GLODAP.(var_type{v})) & GLODAP.(var_type{v}) > 0;
    GLODAP.(var_type{v}) = GLODAP.(var_type{v})(idx);
    GLODAP.([var_type{v} '_grid']) = GLODAP.([var_type{v} '_grid'])(idx);
end
% calculate differences
for v = 1:length(var_type)
    GLODAP.([var_type{v} '_del']) = ...
        GLODAP.(var_type{v}) - GLODAP.([var_type{v} '_grid']);
end

%% show error stats
% pCO2
disp(['Med. delta pCO2 (GLODAP) = ' num2str(round(median(GLODAP.pco2_del,'omitnan'),1))]);
disp(['IQR pCO2 (GLODAP) = ' num2str(round(iqr(GLODAP.pco2_del),1))]);
idx = ~isnan(GLODAP.pco2_grid);
disp(['GLODAP pCO2 ~ RFR-LME pCO2 = ' num2str(round(corr(GLODAP.pco2(idx),GLODAP.pco2_grid(idx)),2))]);
disp(['Propagated pCO2 Err. (GLODAP) = ' num2str(round(mean(GLODAP.pco2_e,'omitnan'),1))]);
disp(['Med. delta pCO2 (CODAP) = ' num2str(round(median(CODAP.pco2_del,'omitnan'),1))]);
disp(['IQR pCO2 (CODAP) = ' num2str(round(iqr(CODAP.pco2_del),1))]);
idx = ~isnan(CODAP.pco2_grid);
disp(['CODAP pCO2 ~ RFR-LME pCO2 = ' num2str(round(corr(CODAP.pco2(idx),CODAP.pco2_grid(idx)),2))]);
disp(['Propagated pCO2 Err. (CODAP) = ' num2str(round(mean(CODAP.pco2_e,'omitnan'),1))]);
disp(' ');
% pH
disp(['Med. delta pH (GLODAP) = ' num2str(round(median(GLODAP.pH_del,'omitnan'),3))]);
disp(['IQR pH (GLODAP) = ' num2str(round(iqr(GLODAP.pH_del),3))]);
idx = ~isnan(GLODAP.pH_grid);
disp(['GLODAP pH ~ RFR-LME pH = ' num2str(round(corr(GLODAP.pH(idx),GLODAP.pH_grid(idx)),2))]);
disp(['Propagated pH Err. (GLODAP) = ' num2str(round(mean(GLODAP.pH_e,'omitnan'),3))]);
disp(['Med. delta pH (CODAP) = ' num2str(round(median(CODAP.pH_del,'omitnan'),3))]);
disp(['IQR pH (CODAP) = ' num2str(round(iqr(CODAP.pH_del),3))]);
idx = ~isnan(CODAP.pH_grid);
disp(['CODAP pH ~ RFR-LME pH = ' num2str(round(corr(CODAP.pH(idx),CODAP.pH_grid(idx)),2))]);
disp(['Propagated pH Err. (CODAP) = ' num2str(round(mean(CODAP.pH_e,'omitnan'),3))]);
disp(' ');
% Omega
disp(['Med. delta OmA (GLODAP) = ' num2str(round(median(GLODAP.OmA_del,'omitnan'),2))]);
disp(['IQR OmA (GLODAP) = ' num2str(round(iqr(GLODAP.OmA_del),2))]);
idx = ~isnan(GLODAP.OmA_grid);
disp(['GLODAP OmA ~ RFR-LME OmA = ' num2str(round(corr(GLODAP.OmA(idx),GLODAP.OmA_grid(idx)),2))]);
disp(['Propagated OmA Err. (GLODAP) = ' num2str(round(mean(GLODAP.OmA_e,'omitnan'),3))]);
disp(['Med. delta OmA (CODAP) = ' num2str(round(median(CODAP.OmA_del,'omitnan'),2))]);
disp(['IQR OmA (CODAP) = ' num2str(round(iqr(CODAP.OmA_del),2))]);
idx = ~isnan(CODAP.OmA_grid);
disp(['CODAP OmA ~ RFR-LME OmA = ' num2str(round(corr(CODAP.OmA(idx),CODAP.OmA_grid(idx)),2))]);
disp(['Propagated OmA Err. (CODAP) = ' num2str(round(mean(CODAP.OmA_e,'omitnan'),3))]);
disp(' ');

%% aggregate deltas on grid
lon_grid = round(min(RFR_LME.lon)):round(max(RFR_LME.lon));
lat_grid = round(min(RFR_LME.lat)):round(max(RFR_LME.lat));
% Determine bin number of each data point
[~,~,Xnum] = histcounts([GLODAP.lon;CODAP.lon],lon_grid);
[~,~,Ynum] = histcounts([GLODAP.lat;CODAP.lat],lat_grid);
subs = [Xnum, Ynum];
sz = [length(lon_grid) length(lat_grid)];
% pCO2
GLODAP_idx = ~isnan(GLODAP.pco2_del);
CODAP_idx = ~isnan(CODAP.pco2_del);
pco2_del_gridded = ...
    accumarray(subs([GLODAP_idx;CODAP_idx],:),...
    [GLODAP.pco2_del(GLODAP_idx);CODAP.pco2_del(CODAP_idx)],...
    sz,@nanmean,NaN);
% pH
GLODAP_idx = ~isnan(GLODAP.pH_del);
CODAP_idx = ~isnan(CODAP.pH_del);
pH_del_gridded = ...
    accumarray(subs([GLODAP_idx;CODAP_idx],:),...
    [GLODAP.pH_del(GLODAP_idx);CODAP.pH_del(CODAP_idx)],...
    sz,@nanmean,NaN);
% omA
GLODAP_idx = ~isnan(GLODAP.OmA_del);
CODAP_idx = ~isnan(CODAP.OmA_del);
OmA_del_gridded = ...
    accumarray(subs([GLODAP_idx;CODAP_idx],:),...
    [GLODAP.OmA_del(GLODAP_idx);CODAP.OmA_del(CODAP_idx)],...
    sz,@nanmean,NaN);

%% plot gridded delta pCO2
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
caxis([-100 100]);
colormap(cmocean('balance','pivot',0));
c.TickLength = 0;
c.Label.String = ['{\itp}CO_{2(discrete)} - {\itp}CO_{2(RFR-LME)}'];
cbarrow;
% plot
pcolorm(double(repmat(lat_grid,length(lon_grid),1)),...
    double(repmat(convert_lon(lon_grid)',1,length(lat_grid))),...
    pco2_del_gridded);
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
exportgraphics(gcf,'Figures/pCO2_delta_eval_map.png');
close

%% plot gridded delta pH
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
caxis([-0.1 0.1]);
colormap(cmocean('balance','pivot',0));
c.TickLength = 0;
c.Label.String = ['pH_{T(discrete)} - pH_{T(RFR-LME)}'];
cbarrow;
% plot
pcolorm(double(repmat(lat_grid,length(lon_grid),1)),...
    double(repmat(convert_lon(lon_grid)',1,length(lat_grid))),...
    pH_del_gridded);
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
exportgraphics(gcf,'Figures/pH_delta_eval_map.png');
close

%% plot gridded delta Omega
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
caxis([-0.2 0.2]);
colormap(cmocean('balance','pivot',0));
c.TickLength = 0;
c.Label.String = ['\Omega_{ar(discrete)} - \Omega_{ar(RFR-LME)}'];
cbarrow;
% plot
pcolorm(double(repmat(lat_grid,length(lon_grid),1)),...
    double(repmat(convert_lon(lon_grid)',1,length(lat_grid))),...
    pH_del_gridded);
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
exportgraphics(gcf,'Figures/OmA_delta_eval_map.png');
%close

%% plot figures
% for v = 2:length(var_type)
%     plot_delta_eval(edges{v},CODAP.(var_type{v}),...
%         CODAP.([var_type{v} '_grid']),var_type{v},var_lab{v},units{v},rounder(v),'CODAP')
% end

%% TA
% CODAP
figure; hold on;
title('TA');
scatter(CODAP.TA,CODAP.TA_grid,'k.');
plot([0 2500],[0 2500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([1500 2500]);
ylim([1500 2500]);
text(1500,2400,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.TA_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.TA_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_TA.png');
% GLODAP
figure; hold on;
title('TA');
scatter(GLODAP.TA,GLODAP.TA_grid,'k.');
plot([1500 2500],[1500 2500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([1500 2500]);
ylim([1500 2500]);
text(1500,2400,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.TA_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.TA_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_TA.png');

%% fCO2
% CODAP
figure; hold on;
title('fCO2');
scatter(CODAP.fco2,CODAP.fco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_fCO2.png');
close
% GLODAP
figure; hold on;
title('fCO2');
scatter(GLODAP.fco2,GLODAP.fco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_fCO2.png');
close
% Both
figure; hold on;
scatter(GLODAP.fco2,GLODAP.fco2_grid,'g.');
scatter(CODAP.fco2,CODAP.fco2_grid,'r.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP/CODAP {\itf}CO_{2}');
ylabel('LME Grid {\itf}CO_{2}');
xlim([0 1200]);
ylim([0 1200]);
text(200,1100,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.fco2_del,[],'omitnan'),1))]);
text(200,900,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_fCO2.png');
close
% Both histogram
figure; hold on;
histogram(GLODAP.fco2_del);
histogram(CODAP.fco2_del);
xlabel('{\itf}CO_{2(discrete)} - {\itf}CO_{2(RFR-LME)}');
ylabel('Counts');
xlim([-500 1000]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_fCO2.png');
close

%% pCO2
% CODAP
figure; hold on;
title('pCO2');
scatter(CODAP.pco2,CODAP.pco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_pCO2.png');
close
% GLODAP
figure; hold on;
title('pCO2');
scatter(GLODAP.pco2,GLODAP.pco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_pCO2.png');
close
% Both
figure; hold on;
scatter(GLODAP.pco2,GLODAP.pco2_grid,'g.');
scatter(CODAP.pco2,CODAP.pco2_grid,'r.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP/CODAP {\itp}CO_{2}');
ylabel('LME Grid {\itp}CO_{2}');
xlim([0 1200]);
ylim([0 1200]);
text(200,1100,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pco2_del,[],'omitnan'),1))]);
text(200,900,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_pCO2.png');
close
% Both histogram
figure; hold on;
histogram(GLODAP.pco2_del);
histogram(CODAP.pco2_del);
plot([0 0],[0 800],'--k','linewidth',2);
xlabel('{\itp}CO_{2(discrete)} - {\itp}CO_{2(RFR-LME)}');
ylabel('Counts');
xlim([-500 1000]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_pCO2.png');
close

%% pH
% CODAP
figure; hold on;
title('pH');
scatter(CODAP.pH,CODAP.pH_grid,'k.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([7.4 8.5]);
ylim([7.4 8.5]);
text(7.5,8.4,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/CODAP_eval_pH.png');
close
% GLODAP
figure; hold on;
title('pH');
scatter(GLODAP.pH,GLODAP.pH_grid,'k.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([7.4 8.6]);
ylim([7.4 8.6]);
text(7.5,8.4,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_eval_pH.png');
close
% Both scatter
figure; hold on;
scatter(GLODAP.pH,GLODAP.pH_grid,'g.');
scatter(CODAP.pH,CODAP.pH_grid,'r.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlim([7.4 8.6]);
ylim([7.4 8.6]);
xlabel('GLODAP/CODAP pH_{T}');
ylabel('LME Grid pH_{T}');
text(7.5,8.5,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pH_del,[],'omitnan'),3))]);
text(7.5,8.42,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_pH.png');
% Both histogram
figure; hold on;
histogram(GLODAP.pH_del);
histogram(CODAP.pH_del);
plot([0 0],[0 900],'--k','linewidth',2);
xlabel('pH_{T(discrete)} - pH_{T(RFR-LME)}');
ylabel('Counts');
xlim([-1 1]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_pH.png');
close

%% OmA
% CODAP
figure; hold on;
title('\Omega_{A}');
scatter(CODAP.OmA,CODAP.OmA_grid,'k.');
plot([0 5],[0 5],'k-');
xlabel('CODAP');
ylabel('LME Grid');
text(0.5,4.5,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/CODAP_eval_OmA.png');
close
% GLODAP
figure; hold on;
title('\Omega_{A}');
scatter(GLODAP.OmA,GLODAP.OmA_grid,'k.');
plot([0 5],[0 5],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
text(0.5,4.5,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_eval_OmA.png');
close
% Both scatter
figure; hold on;
scatter(GLODAP.OmA,GLODAP.OmA_grid,'g.');
scatter(CODAP.OmA,CODAP.OmA_grid,'r.');
plot([0 5],[0 5],'k-');
ylim([0 5]); xlim([0 5]);
xlabel('GLODAP/CODAP \Omega_{A}');
ylabel('LME Grid \Omega_{A}');
text(0.5,4.7,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.OmA_del,[],'omitnan'),3))]);
text(0.5,4.4,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_OmA.png');
% Both histogram
figure; hold on;
histogram(GLODAP.OmA_del);
histogram(CODAP.OmA_del);
plot([0 0],[0 450],'--k','linewidth',2);
xlabel('\Omega_{ar(discrete)} - \Omega_{ar(RFR-LME)}');
ylabel('Counts');
xlim([-2.5 2.5]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_OmA.png');
close

% clean up
clear
close all