% Convert MATLAB files to NetCDFs
% 
% This script takes regional MATLAB files and converts them to full region
% NetCDF files.
% 
% Written by J.D. Sharp: 1/30/23
% Last updated by J.D. Sharp: 1/30/23
% 

%% this script defines the bounds of the eighteen LMEs
define_regions

%% full data grid structure
load('Data/socat_gridded','SOCAT_grid');
US_LME_RFR.lim = SOCAT_grid.lim;
US_LME_RFR.dim = SOCAT_grid.dim;
US_LME_RFR.lon = SOCAT_grid.lon;
US_LME_RFR.lat = SOCAT_grid.lat;
US_LME_RFR.month = SOCAT_grid.month;
US_LME_RFR.year = SOCAT_grid.year;
clear SOCAT_grid

%% Pre-allocate full grid of predictors
Vars_pred = {'SSS' 'SSH' 'SST' 'IceC' 'CHL' 'WindSpeed' 'Bathy' ... 
    'MLD' 'mslp' 'pCO2_atm'};
for v = 1:length(Vars_pred)
   if ~strcmp(Vars_pred{v},'Bathy')
       US_LME_RFR.(Vars_pred{v}) = ...
           nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
   else
       US_LME_RFR.(Vars_pred{v}) = ...
           nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y);
   end
end

%% Pre-allocate full grid of OA Indicators
Vars_OA = {'fCO2' 'DIC' 'TA' 'uTA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
for v = 1:length(Vars_OA)
   US_LME_RFR.(Vars_OA{v}) = ...
       nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
end

%% Add regional variables to full grid
for n = 1:length(region)

    %% load gridded fCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');

    %% Add regional predictors to full grid
    idx_lon = US_LME_RFR.lon >= min(Preds_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(Preds_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(Preds_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(Preds_grid.(region{n}).lat);
    for v = 1:length(Vars_pred)
        tmp_var = US_LME_RFR.(Vars_pred{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            Preds_grid.(region{n}).(Vars_pred{v});
        idx = ~isnan(tmp_var);
        US_LME_RFR.(Vars_pred{v})(idx) = tmp_var(idx);
    end

    %% Add OA indicators to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);
    for v = 1:length(Vars_OA)
        tmp_var = US_LME_RFR.(Vars_OA{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            OAI_grid.(region{n}).(Vars_OA{v});
        idx = ~isnan(tmp_var);
        US_LME_RFR.(Vars_OA{v})(idx) = tmp_var(idx);
    end

    %% clean up
    clear idx_lon idx_lat idx v Preds_grid OAI_grid tmp_var

end

%% save predictors as NetCDF
filename = ['Data/US_LME_RFR_Preds_' date '.nc'];
time = datenum([repmat(1998,US_LME_RFR.dim.z,1) US_LME_RFR.month, ...
    zeros(US_LME_RFR.dim.z,1)]);
nccreate(filename,'Lon','Dimensions',{'Lon',US_LME_RFR.dim.x});
ncwrite(filename,'Lon',US_LME_RFR.lon);
ncwriteatt(filename,'Lon','_CoordinateAxisType','Lon');
nccreate(filename,'Lat','Dimensions',{'Lat',US_LME_RFR.dim.y});
ncwrite(filename,'Lat',US_LME_RFR.lat);
ncwriteatt(filename,'Lat','_CoordinateAxisType','Lat');
nccreate(filename,'Time','Dimensions',{'Time',US_LME_RFR.dim.z});
ncwrite(filename,'Time',time);
ncwriteatt(filename,'Time','_CoordinateAxisType','Time');
for v = 1:length(Vars_pred)
    if ~strcmp(Vars_pred{v},'Bathy')
        nccreate(filename,Vars_pred{v},'Dimensions',{'Lon',US_LME_RFR.dim.x,'Lat',US_LME_RFR.dim.y,'Time',US_LME_RFR.dim.z});
    else
        nccreate(filename,Vars_pred{v},'Dimensions',{'Lon',US_LME_RFR.dim.x,'Lat',US_LME_RFR.dim.y});
    end
    ncwrite(filename,Vars_pred{v},US_LME_RFR.(Vars_pred{v}));
end

%% save OA indicators as NetCDF
filename = ['Data/US_LME_RFR_Inds_' date '.nc'];
time = datenum([repmat(1998,US_LME_RFR.dim.z,1) US_LME_RFR.month, ...
    zeros(US_LME_RFR.dim.z,1)]);
nccreate(filename,'Lon','Dimensions',{'Lon',US_LME_RFR.dim.x});
ncwrite(filename,'Lon',US_LME_RFR.lon);
ncwriteatt(filename,'Lon','_CoordinateAxisType','Lon');
nccreate(filename,'Lat','Dimensions',{'Lat',US_LME_RFR.dim.y});
ncwrite(filename,'Lat',US_LME_RFR.lat);
ncwriteatt(filename,'Lat','_CoordinateAxisType','Lat');
nccreate(filename,'Time','Dimensions',{'Time',US_LME_RFR.dim.z});
ncwrite(filename,'Time',time);
ncwriteatt(filename,'Time','_CoordinateAxisType','Time');
for v = 1:length(Vars_OA)
    nccreate(filename,Vars_OA{v},'Dimensions',{'Lon',US_LME_RFR.dim.x,'Lat',US_LME_RFR.dim.y,'Time',US_LME_RFR.dim.z});
    ncwrite(filename,Vars_OA{v},US_LME_RFR.(Vars_OA{v}));
end

%% plot for sanity
figure;
worldmap([US_LME_RFR.lim.latmin US_LME_RFR.lim.latmax],...
    [US_LME_RFR.lim.lonmin US_LME_RFR.lim.lonmax]);
pcolorm(US_LME_RFR.lat,US_LME_RFR.lon,mean(US_LME_RFR.fCO2,3,'omitnan')');
close

%% Clean up
clear US_LME_RFR Vars_OA Vars_pred v n time filename
