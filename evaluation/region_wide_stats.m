%% load US LME data
date = '02-Aug-2023';
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

%% Load SOCAT grid
load('Data/socat_gridded_2023','SOCAT_grid');

%% calculate area-weighted time series
area_weights = SOCAT_grid.area_km2.*SOCAT_grid.percent_sea;
fCO2_dom_mean = nan(SOCAT_grid.dim.z,1);
pH_dom_mean = nan(SOCAT_grid.dim.z,1);
OmA_dom_mean = nan(SOCAT_grid.dim.z,1);
for t = 1:SOCAT_grid.dim.z
        % remove ice-filled cells from area weights
        area_weights(isnan(LME_RFR.fco2(:,:,t))) = NaN;
        % calculate area-weighted means
        pCO2_dom_mean(t) = ...
            squeeze(sum(sum(LME_RFR.pco2(:,:,t).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        pH_dom_mean(t) = ...
            squeeze(sum(sum(LME_RFR.pH(:,:,t).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        OmA_dom_mean(t) = ...
            squeeze(sum(sum(LME_RFR.OmA(:,:,t).*...
                area_weights,1,'omitnan'),2,'omitnan'))./...
                squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
end

%% calculate trend
% pCO2
[yf,yr,x,err] = ...
    leastsq2(LME_RFR.time,pCO2_dom_mean,LME_RFR.time(1),2,[365/2 365]);
x(2)*365
[acov,acor,lag,dof] = autocov2(LME_RFR.time,yr,365*3);
edof = dof - 2; % subtract number of parameters to get effective dof
tr_uncer = ... % scale uncertainty using edof
    err(2)*(sqrt(length(LME_RFR.time./sqrt(edof))))*365
% pH
[yf,yr,x,err] = ...
    leastsq2(LME_RFR.time,pH_dom_mean,LME_RFR.time(1),2,[365/2 365]);
x(2)*365
[acov,acor,lag,dof] = autocov2(LME_RFR.time,yr,365*3);
edof = dof - 2; % subtract number of parameters to get effective dof
tr_uncer = ... % scale uncertainty using edof
    err(2)*(sqrt(length(LME_RFR.time./sqrt(edof))))*365
% OmA
[yf,yr,x,err] = ...
    leastsq2(LME_RFR.time,OmA_dom_mean,LME_RFR.time(1),2,[365/2 365]);
x(2)*365
[acov,acor,lag,dof] = autocov2(LME_RFR.time,yr,365*3);
edof = dof - 2; % subtract number of parameters to get effective dof
tr_uncer = ... % scale uncertainty using edof
    err(2)*(sqrt(length(LME_RFR.time./sqrt(edof))))*365
