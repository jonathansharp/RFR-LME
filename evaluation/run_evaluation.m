% Run scripts to evaluate US LME OA indicators

%% load CODAP data
load('CODAP_NA/CODAP_NA_v2020_G2format.mat')
idx = C1.pressure <= 20 & C1.tco2f == 2 & C1.talkf == 2 & ...
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
CODAP.pH = C1.phtsinsitutp(idx);
CODAP.OmA = C1.aragonite(idx);
CODAP.OmC = C1.calcite(idx);
CODAP.RF = C1.revelle(idx);
clear C1 idx Headers Units

%% calculate fCO2, H, and CO3
carb = CO2SYS(CODAP.TA,CODAP.DIC,1,2,CODAP.sal,CODAP.tmp,NaN,CODAP.prs,...
    NaN,CODAP.sil,CODAP.phos,0,0,1,10,1,2,2);
CODAP.fco2 = carb(:,5);
CODAP.H = carb(:,15).*10^3;
CODAP.CO3 = carb(:,7);
clear carb

%% load US LME data
LME_RFR.lat = ncread('Data/US_LME_RFR_Inds.nc','Lat');
LME_RFR.lon = ncread('Data/US_LME_RFR_Inds.nc','Lon');
LME_RFR.time = ncread('Data/US_LME_RFR_Inds.nc','Time');
LME_RFR.fco2 = ncread('Data/US_LME_RFR_Inds.nc','fCO2');
LME_RFR.TA = ncread('Data/US_LME_RFR_Inds.nc','TA');
% LME_RFR.DIC = ncread('Data/US_LME_RFR_Inds.nc','DIC');
LME_RFR.pH = ncread('Data/US_LME_RFR_Inds.nc','pH');
LME_RFR.OmA = ncread('Data/US_LME_RFR_Inds.nc','OmA');
LME_RFR.OmC = ncread('Data/US_LME_RFR_Inds.nc','OmC');
LME_RFR.H = ncread('Data/US_LME_RFR_Inds.nc','H');
LME_RFR.CO3 = ncread('Data/US_LME_RFR_Inds.nc','CO3');
LME_RFR.RF = ncread('Data/US_LME_RFR_Inds.nc','RF');

%% co-locate each CODAP point with an LME-RFR grid cell
% variable information
edges = {1900:10:2300;300:5:500;2000:10:2400;7.8:0.0125:8.3;...
    1:0.1:5;1:0.1:5;7:0.2:14;100:2.5:250;9:0.1:17};
var_type = {'DIC' 'fco2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 3 2 2 1 1 2];
% pre-allocate
for v = 1:length(var_type)
    CODAP.([var_type{v} '_grid']) = nan(size(CODAP.(var_type{v})));
end
% identify and fill gridded values
for n = 1:length(CODAP.lat)
    idx_lon = find(abs(LME_RFR.lon - CODAP.lon(n)) == ...
        min(abs(LME_RFR.lon - CODAP.lon(n))));
    idx_lat = find(abs(LME_RFR.lat - CODAP.lat(n)) == ...
        min(abs(LME_RFR.lat - CODAP.lat(n))));
    idx_time = find(abs(LME_RFR.time - CODAP.time(n)) == ...
        min(abs(LME_RFR.time - CODAP.time(n))));
    for v = 2:length(var_type)
        CODAP.([var_type{v} '_grid'])(n) = ...
            LME_RFR.(var_type{v})(idx_lon(1),idx_lat(1),idx_time(1));
    end
end

%% plot figures
for v = 2:length(var_type)
    plot_delta_eval(edges{v},CODAP.(var_type{v}),...
        CODAP.([var_type{v} '_grid']),var_type{v},var_lab{v},units{v},rounder(v),'CODAP')
end

%% load GLODAP data
load('GLODAPv2.2022/GLODAPv2.2022_Merged_Master_File.mat')
GLODAP.lat = G2latitude;
GLODAP.lon = G2longitude;
GLODAP.idx = G2talkf == 2 & G2tco2f == 2;
clearvars -except GLODAP