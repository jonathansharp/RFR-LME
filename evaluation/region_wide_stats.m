%% load US LME data
% set date
date = '08-May-2024';
new = 1;
% load indicators
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
RFR_LME.pco2_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pCO2.nc'],'u_pco2');
RFR_LME.fco2_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_fCO2.nc'],'u_fco2');
RFR_LME.TA_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_TA.nc'],'u_ta');
RFR_LME.DIC_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_DIC.nc'],'u_dic');
RFR_LME.pH_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_pH.nc'],'u_ph');
RFR_LME.OmA_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_OmA.nc'],'u_om_a');
RFR_LME.OmC_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_OmC.nc'],'u_om_c');
RFR_LME.H_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_H.nc'],'u_h');
RFR_LME.CO3_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_CO3.nc'],'u_co3');
RFR_LME.RF_e = ncread(['Data/NetCDFs_' date '/US_RFR_LME_RF.nc'],'u_rf');
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

%% Load SOCAT grid
load('Data/socat_gridded_2022','SOCAT_grid');

%% define variables
var_type_1 = {'DIC' 'pCO2' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_type_2 = {'DIC' 'pco2' 'fco2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
dom_mean = struct();
stats_table = nan(6,length(var_type_1));

%% calculate area-weighted time series

% establish area_weights
area_weights = SOCAT_grid.area_km2.*SOCAT_grid.percent_sea;
area_weights = repmat(area_weights,1,1,12);
% establish climatological index of only cells that are always
spatial_index = nan(SOCAT_grid.dim.x,SOCAT_grid.dim.y,12);
for m = 1:12
    spatial_index(:,:,m) = ...
        sum(~isnan(RFR_LME.fco2(:,:,m:12:end)),3) == ...
            SOCAT_grid.dim.z/12;
end
area_weights(~spatial_index) = NaN;
area_weights = repmat(area_weights,1,1,SOCAT_grid.dim.z/12);

% calculate area-weighted means
for v = 1:length(var_type_1)
        dom_mean.(var_type_1{v}) = ...
            squeeze(sum(sum(RFR_LME.(var_type_2{v}).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        dom_mean.([var_type_1{v} '_e']) = ...
            squeeze(sum(sum(RFR_LME.([var_type_2{v} '_e']).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
end

%% calculate stats
for v = 1:length(var_type_1)
    stats_table(1,v) = mean(dom_mean.(var_type_1{v})); % long-term mean
    stats_table(2,v) = mean(dom_mean.([var_type_1{v} '_e'])); % long-term mean error
    % determine climatology
    clim = nan(12,1);
    for m = 1:12
        clim(m) = mean(dom_mean.(var_type_1{v})(m:12:end),'omitnan');
    end
    % calculate amplitude
    if sum(~isnan(clim)) > 8
        amp = max(clim) - min(clim);
    else
        amp = NaN;
    end
    stats_table(3,v) = amp; % amplitude
    [yf,yr,x,err] = ...
        leastsq2(RFR_LME.time,dom_mean.(var_type_1{v}),RFR_LME.time(1),3,[365/3 365/2 365]);
    stats_table(4,v) = std(yr); % iav
    stats_table(5,v) = x(2)*365; % trend
    [acov,acor,lag,dof] = autocov2(RFR_LME.time,yr,365*3);
    edof = dof - 8; % subtract number of parameters to get effective dof
    %edof<1 = 1;
    stats_table(6,v) = ... % scale uncertainty using edof
        err(2)*(sqrt(length(RFR_LME.time./sqrt(edof))))*365;
end

%% save table
stats_table = array2table(stats_table,'RowNames',...
    {'Mean' 'Mean_e' 'Amp' 'IAV' 'Trend' 'Trend Uncer.'},'VariableNames',var_type_1);
writetable(stats_table,['IndsAndStats/GlobalStatsTable-' date '.xls'],'WriteRowNames',true);
clear
