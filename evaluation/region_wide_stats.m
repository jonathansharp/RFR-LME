%% load US LME data
% set date
date = '11-Jan-2024';
new = 1;
% load indicators
if new == 1
% New NetCDF files
LME_RFR.lat = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pCO2.nc'],'lat');
LME_RFR.lon = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pCO2.nc'],'lon');
LME_RFR.time = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pCO2.nc'],'time');
LME_RFR.pco2 = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pCO2.nc'],'pco2');
LME_RFR.fco2 = ncread(['Data/NetCDFs_' date '/US_LME_RFR_fCO2.nc'],'fco2');
LME_RFR.TA = ncread(['Data/NetCDFs_' date '/US_LME_RFR_TA.nc'],'ta');
LME_RFR.DIC = ncread(['Data/NetCDFs_' date '/US_LME_RFR_DIC.nc'],'dic');
LME_RFR.pH = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pH.nc'],'ph');
LME_RFR.OmA = ncread(['Data/NetCDFs_' date '/US_LME_RFR_OmA.nc'],'om_a');
LME_RFR.OmC = ncread(['Data/NetCDFs_' date '/US_LME_RFR_OmC.nc'],'om_c');
LME_RFR.H = ncread(['Data/NetCDFs_' date '/US_LME_RFR_H.nc'],'h');
LME_RFR.CO3 = ncread(['Data/NetCDFs_' date '/US_LME_RFR_CO3.nc'],'co3');
LME_RFR.RF = ncread(['Data/NetCDFs_' date '/US_LME_RFR_RF.nc'],'rf');
LME_RFR.pco2_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pCO2.nc'],'u_pco2');
LME_RFR.fco2_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_fCO2.nc'],'u_fco2');
LME_RFR.TA_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_TA.nc'],'u_ta');
LME_RFR.DIC_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_DIC.nc'],'u_dic');
LME_RFR.pH_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_pH.nc'],'u_ph');
LME_RFR.OmA_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_OmA.nc'],'u_om_a');
LME_RFR.OmC_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_OmC.nc'],'u_om_c');
LME_RFR.H_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_H.nc'],'u_h');
LME_RFR.CO3_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_CO3.nc'],'u_co3');
LME_RFR.RF_e = ncread(['Data/NetCDFs_' date '/US_LME_RFR_RF.nc'],'u_rf');
else
% Old NetCDF files
LME_RFR.lat = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Lat');
LME_RFR.lon = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Lon');
LME_RFR.time = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Time');
LME_RFR.pco2 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'pCO2');
LME_RFR.fco2 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'fCO2');
LME_RFR.TA = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'TA');
LME_RFR.DIC = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'DIC');
LME_RFR.pH = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'pH');
LME_RFR.OmA = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'OmA');
LME_RFR.OmC = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'OmC');
LME_RFR.H = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'H');
LME_RFR.CO3 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'CO3');
LME_RFR.RF = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'RF');
end

%% Load SOCAT grid
load('Data/socat_gridded_2023','SOCAT_grid');

%% define variables
var_type_1 = {'DIC' 'pCO2' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_type_2 = {'DIC' 'pco2' 'fco2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
dom_mean = struct();
stats_table = nan(6,length(var_type_1));

%% calculate area-weighted time series
% define area weights
area_weights = SOCAT_grid.area_km2.*SOCAT_grid.percent_sea;
area_weights = repmat(area_weights,1,1,SOCAT_grid.dim.z);
% remove ice-filled cells from area weights
area_weights(isnan(LME_RFR.fco2)) = NaN;
for v = 1:length(var_type_1)
        % calculate area-weighted means
        dom_mean.(var_type_1{v}) = ...
            squeeze(sum(sum(LME_RFR.(var_type_2{v}).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        dom_mean.([var_type_1{v} '_e']) = ...
            squeeze(sum(sum(LME_RFR.([var_type_2{v} '_e']).*...
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
        leastsq2(LME_RFR.time,dom_mean.(var_type_1{v}),LME_RFR.time(1),3,[365/3 365/2 365]);
    stats_table(4,v) = std(yr); % iav
    stats_table(5,v) = x(2)*365; % trend
    [acov,acor,lag,dof] = autocov2(LME_RFR.time,yr,365*3);
    edof = dof - 8; % subtract number of parameters to get effective dof
    %edof<1 = 1;
    stats_table(6,v) = ... % scale uncertainty using edof
        err(2)*(sqrt(length(LME_RFR.time./sqrt(edof))))*365;
end

%% save table
stats_table = array2table(stats_table,'RowNames',...
    {'Mean' 'Mean_e' 'Amp' 'IAV' 'Trend' 'Trend Uncer.'},'VariableNames',var_type_1);
writetable(stats_table,['IndsAndStats/GlobalStatsTable-' date '.xls'],'WriteRowNames',true);
clear
