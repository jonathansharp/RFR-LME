% Convert MATLAB files to NetCDFs
% 
% This script takes regional MATLAB files and converts them to full region
% NetCDF files.
% 
% Written by J.D. Sharp: 1/30/23
% Last updated by J.D. Sharp: 1/17/24
% 

%% this script defines the bounds of the eleven LMEs
define_regions_eiwg

%% full data grid structure
load('Data/socat_gridded','SOCAT_grid');
US_RFR_LME.lim = SOCAT_grid.lim;
US_RFR_LME.dim = SOCAT_grid.dim;
US_RFR_LME.lon = single(SOCAT_grid.lon);
US_RFR_LME.lat = single(SOCAT_grid.lat);
US_RFR_LME.month = single(SOCAT_grid.month);
US_RFR_LME.year = single(SOCAT_grid.year);
clear SOCAT_grid

%% Pre-allocate full grid of LME numbers
US_RFR_LME.LME = ...
    single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y));

%% Pre-allocate full grid of predictors
Vars_pred = {'SSS' 'SSH' 'SST' 'IceC' 'CHL' 'WindSpeed' 'Bathy' ... 
    'MLD' 'mslp' 'pCO2_atm'};
varnames_pred = {'sss' 'ssh' 'sst' 'ice_c' 'chl' 'ws' 'bottom_depth' ... 
    'mld' 'p' 'pco2_a'};
names_pred = {'Sea Surface Salinity' 'Sea Surface Height' 'Sea Surface Temperature' ...
    'Ice Concentration' 'Chlorophyll-a Concentration' 'Wind Speed' 'Bathymetry' ...
    'Mixed Layer Depth' 'Sea Level Pressure' 'Atmosphreric pCO2'};
units_pred = {'' 'millimeters' 'degrees Celcius' 'percent of area' ...
    'milligrams per meter squared' 'meters per second' 'meters' 'meters' ...
    '' 'microatmospheres'};
for v = 1:length(Vars_pred)
   if ~strcmp(Vars_pred{v},'Bathy') && ~strcmp(Vars_pred{v},'Dist')
       US_RFR_LME.(Vars_pred{v}) = ...
           single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y,US_RFR_LME.dim.z));
   else
       US_RFR_LME.(Vars_pred{v}) = ...
           single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y));
   end
end

%% Pre-allocate full grid of OA Indicators and Uncertainties
% variable names for matlab structure
Vars_OA = {'fCO2' 'pCO2' 'DIC' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
uVars_OA = {'ufCO2' 'upCO2' 'uDIC' 'uTA' 'upH' 'uOmA' 'uOmC' 'uH' 'uCO3' 'uRF'};
% variable names for NetCDF file
varnames_OA = {'fco2' 'pco2' 'dic' 'ta' 'ph' 'om_a' 'om_c' 'h' 'co3' 'rf'};
varnames_uOA = {'u_fco2' 'u_pco2' 'u_dic' 'u_ta' 'u_ph' 'u_om_a' 'u_om_c' ...
    'u_h' 'u_co3' 'u_rf'};
% long variable names
names_OA = {'Sea Surface CO2 Fugacity' 'Sea Surface CO2 Partial Pressure' ...
    'Sea Surface Dissolved Inorganic Carbon' 'Sea Surface Total Alkalinity' ...
    'Sea Surface pH (total scale)' 'Sea Surface Aragonite Saturation State' ...
    'Sea Surface Calcite Saturation State' 'Sea Surface Hydrogen Ion Concentration' ...
    'Sea Surface Carbonate Ion Concentration' 'Sea Surface Revelle Factor'};
% variable units
units_OA = {'microatmospheres' 'microatmospheres' 'micromoles per kilogram' ...
    'micromoles per kilogram' '' '' '' 'nanomoles per kilogram' 'micromoles per kilogram' ''};
% pre-allocate variables grid
for v = 1:length(Vars_OA)
   US_RFR_LME.(Vars_OA{v}) = ...
       single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y,US_RFR_LME.dim.z));
end
% pre-allocate uncertainties grid
for v = 1:length(uVars_OA)
   US_RFR_LME.(uVars_OA{v}) = ...
       single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y,US_RFR_LME.dim.z));
end

%% Add regional variables to full grid
for n = 1:length(region)

    %% load gridded fCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');

    %% define LME mask
    idx_lon = US_RFR_LME.lon >= min(Preds_grid.(region{n}).lon) & ...
        US_RFR_LME.lon <= max(Preds_grid.(region{n}).lon);
    idx_lat = US_RFR_LME.lat >= min(Preds_grid.(region{n}).lat) & ...
        US_RFR_LME.lat <= max(Preds_grid.(region{n}).lat);

    %% Add regional predictors to full grid
    for v = 1:length(Vars_pred)
        tmp_var = US_RFR_LME.(Vars_pred{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            Preds_grid.(region{n}).(Vars_pred{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(Vars_pred{v})(idx) = tmp_var(idx);
    end    

    %% Add OA indicators to full grid
    for v = 1:length(Vars_OA)
        tmp_var = US_RFR_LME.(Vars_OA{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            OAI_grid.(region{n}).(Vars_OA{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(Vars_OA{v})(idx) = tmp_var(idx);
    end

    %% Add OA indicator uncertainties to full grid
    for v = 1:length(uVars_OA)
        tmp_var = US_RFR_LME.(uVars_OA{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            OAI_grid.(region{n}).(uVars_OA{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(uVars_OA{v})(idx) = tmp_var(idx);
    end

    %% add LME number to grid
    tmp_idx = false(size(US_RFR_LME.LME));
    tmp_idx(idx_lon,idx_lat) = ...
        any(OAI_grid.(region{n}).idxspc,3);
    tmp_var = US_RFR_LME.LME;
    tmp_var(tmp_idx) = n;
    idx = ~isnan(tmp_var);
    US_RFR_LME.LME(idx) = tmp_var(idx);

    %% clean up
    clear idx_lon idx_lat idx v Preds_grid OAI_grid tmp_var

end

%% save predictors as NetCDFs
% process time
time = datenum([repmat(1998,US_RFR_LME.dim.z,1),...
    double(US_RFR_LME.month+0.5),repmat(15,US_RFR_LME.dim.z,1)]);
% create folder
if ~isfolder([pwd '/Data/NetCDFs_' date]); mkdir(['Data/NetCDFs_' date]); end
for v = 1:length(Vars_pred)
    % define file name
    filename = ['Data/NetCDFs_' date '/US_RFR_LME_' Vars_pred{v} '.nc'];
    % delete old file if it exists
    if isfile(filename); delete(filename); end
    % longitude
    nccreate(filename,'lon','Dimensions',{'lon',US_RFR_LME.dim.x},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lon',US_RFR_LME.lon);
    %ncwriteatt(filename,'lon','bounds','lon_bounds');
    ncwriteatt(filename,'lon','units','degrees_east');
    ncwriteatt(filename,'lon','axis','X');
    ncwriteatt(filename,'lon','long_name','longitude');
    ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
    % latitude
    nccreate(filename,'lat','Dimensions',{'lat',US_RFR_LME.dim.y},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lat',US_RFR_LME.lat);
    %ncwriteatt(filename,'lat','bounds','lat_bounds');
    ncwriteatt(filename,'lat','units','degrees_north');
    ncwriteatt(filename,'lat','axis','Y');
    ncwriteatt(filename,'lat','long_name','latitude');
    ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
    % time
    nccreate(filename,'time','Dimensions',{'time',US_RFR_LME.dim.z},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'time',time);
    %ncwriteatt(filename,'time','bounds','time_bounds');
    ncwriteatt(filename,'time','units','days since 1-Jan-0000');
    ncwriteatt(filename,'time','axis','T');
    ncwriteatt(filename,'time','long_name','time');
    ncwriteatt(filename,'time','_CoordinateAxisType','Time');
    % variable
    if ~strcmp(Vars_pred{v},'Bathy') && ~strcmp(Vars_pred{v},'Dist')
        nccreate(filename,varnames_pred{v},'Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
            'time',US_RFR_LME.dim.z},'DataType','single','FillValue',NaN);
    else
        nccreate(filename,varnames_pred{v},'Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y},...
            'DataType','single','FillValue',NaN);
    end
    ncwrite(filename,varnames_pred{v},US_RFR_LME.(Vars_pred{v}));
    ncwriteatt(filename,varnames_pred{v},'units',units_pred{v});
    ncwriteatt(filename,varnames_pred{v},'long_name',names_pred{v});
    % LME number
    nccreate(filename,'lme_num','Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y},...
            'DataType','single','FillValue',NaN);
    ncwrite(filename,'lme_num',US_RFR_LME.LME);
    ncwriteatt(filename,'lme_num','reference_table','lme_ref');
    ncwriteatt(filename,'lme_num','long_name','Large Marine Ecosystem Number');
    ncwriteatt(filename,'lme_num','information','See "lme_ref" for the order of LMEs associated with each LME number.');
    % LME reference
    nccreate(filename,'lme_ref','Dimensions',{'a',11,'b',4},'DataType','char');
    ncwrite(filename,'lme_ref',char(region));
    ncwriteatt(filename,'lme_num','long_name','Large Marine Ecosystem Reference Table');
    ncwriteatt(filename,'lme_num','information','This variable provides the order of LME abbreviations associated with each LME number.');
end

%% save OA indicators and uncertainties as NetCDFs
% process time
time = datenum([repmat(1998,US_RFR_LME.dim.z,1),...
    double(US_RFR_LME.month+0.5),repmat(15,US_RFR_LME.dim.z,1)]);
% create folder
if ~isfolder([pwd '/Data/NetCDFs_' date]); mkdir(['Data/NetCDFs_' date]); end
for v = 1:length(Vars_OA)
    % define file name
    filename = ['Data/NetCDFs_' date '/US_RFR_LME_' Vars_OA{v} '.nc'];
    % delete old file if it exists
    if isfile(filename); delete(filename); end
    % longitude
    nccreate(filename,'lon','Dimensions',{'lon',US_RFR_LME.dim.x},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lon',US_RFR_LME.lon);
    %ncwriteatt(filename,'lon','bounds','lon_bounds');
    ncwriteatt(filename,'lon','units','degrees_east');
    ncwriteatt(filename,'lon','axis','X');
    ncwriteatt(filename,'lon','long_name','longitude');
    ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
    % latitude
    nccreate(filename,'lat','Dimensions',{'lat',US_RFR_LME.dim.y},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lat',US_RFR_LME.lat);
    %ncwriteatt(filename,'lat','bounds','lat_bounds');
    ncwriteatt(filename,'lat','units','degrees_north');
    ncwriteatt(filename,'lat','axis','Y');
    ncwriteatt(filename,'lat','long_name','latitude');
    ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
    % time
    nccreate(filename,'time','Dimensions',{'time',US_RFR_LME.dim.z},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'time',time);
    %ncwriteatt(filename,'time','bounds','time_bounds');
    ncwriteatt(filename,'time','units','days since 1-Jan-0000');
    ncwriteatt(filename,'time','axis','T');
    ncwriteatt(filename,'time','long_name','time');
    ncwriteatt(filename,'time','_CoordinateAxisType','Time');
    % variable
    nccreate(filename,varnames_OA{v},'Dimensions',...
        {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
        'time',US_RFR_LME.dim.z},'DataType','single','FillValue',NaN);
    ncwrite(filename,varnames_OA{v},US_RFR_LME.(Vars_OA{v}));
    ncwriteatt(filename,varnames_OA{v},'units',units_OA{v});
    ncwriteatt(filename,varnames_OA{v},'long_name',names_OA{v});
    % variable uncertainty
    nccreate(filename,varnames_uOA{v},'Dimensions',...
        {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
        'time',US_RFR_LME.dim.z},'DataType','single','FillValue',NaN);
    ncwrite(filename,varnames_uOA{v},US_RFR_LME.(uVars_OA{v}));
    ncwriteatt(filename,varnames_uOA{v},'units',units_OA{v});
    ncwriteatt(filename,varnames_uOA{v},'long_name',['Uncertainty in ' names_OA{v}]);
    % LME number
    nccreate(filename,'lme_num','Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y},...
            'DataType','single','FillValue',NaN);
    ncwrite(filename,'lme_num',US_RFR_LME.LME);
    ncwriteatt(filename,'lme_num','reference_table','lme_ref');
    ncwriteatt(filename,'lme_num','long_name','Large Marine Ecosystem Number');
    ncwriteatt(filename,'lme_num','information','See "lme_ref" for the order of LMEs associated with each LME number.');
    % LME reference
    nccreate(filename,'lme_ref','Dimensions',{'a',11,'b',4},'DataType','char');
    ncwrite(filename,'lme_ref',char(region));
    ncwriteatt(filename,'lme_num','long_name','Large Marine Ecosystem Reference Table');
    ncwriteatt(filename,'lme_num','information','This variable provides the order of LME abbreviations associated with each LME number.');
end

% %% save predictors as NetCDFs
% filename = ['Data/US_RFR_LME_Preds_' date '.nc'];
% time = datenum([repmat(1998,US_RFR_LME.dim.z,1) US_RFR_LME.month+0.5, ...
%     repmat(15,US_RFR_LME.dim.z,1)]);
% nccreate(filename,'Lon','Dimensions',{'Lon',US_RFR_LME.dim.x});
% ncwrite(filename,'Lon',US_RFR_LME.lon);
% ncwriteatt(filename,'Lon','_CoordinateAxisType','Lon');
% nccreate(filename,'Lat','Dimensions',{'Lat',US_RFR_LME.dim.y});
% ncwrite(filename,'Lat',US_RFR_LME.lat);
% ncwriteatt(filename,'Lat','_CoordinateAxisType','Lat');
% nccreate(filename,'Time','Dimensions',{'Time',US_RFR_LME.dim.z});
% ncwrite(filename,'Time',time);
% ncwriteatt(filename,'Time','_CoordinateAxisType','Time');
% for v = 1:length(Vars_pred)
%     if ~strcmp(Vars_pred{v},'Bathy') && ~strcmp(Vars_pred{v},'Dist')
%         nccreate(filename,Vars_pred{v},'Dimensions',{'Lon',US_RFR_LME.dim.x,'Lat',US_RFR_LME.dim.y,'Time',US_RFR_LME.dim.z});
%     else
%         nccreate(filename,Vars_pred{v},'Dimensions',{'Lon',US_RFR_LME.dim.x,'Lat',US_RFR_LME.dim.y});
%     end
%     ncwrite(filename,Vars_pred{v},US_RFR_LME.(Vars_pred{v}));
% end
% 
% %% save OA indicators as NetCDFs
% filename = ['Data/US_RFR_LME_Inds_' date '.nc'];
% time = datenum([repmat(1998,US_RFR_LME.dim.z,1) US_RFR_LME.month+0.5, ...
%     repmat(15,US_RFR_LME.dim.z,1)]);
% nccreate(filename,'Lon','Dimensions',{'Lon',US_RFR_LME.dim.x});
% ncwrite(filename,'Lon',US_RFR_LME.lon);
% ncwriteatt(filename,'Lon','_CoordinateAxisType','Lon');
% nccreate(filename,'Lat','Dimensions',{'Lat',US_RFR_LME.dim.y});
% ncwrite(filename,'Lat',US_RFR_LME.lat);
% ncwriteatt(filename,'Lat','_CoordinateAxisType','Lat');
% nccreate(filename,'Time','Dimensions',{'Time',US_RFR_LME.dim.z});
% ncwrite(filename,'Time',time);
% ncwriteatt(filename,'Time','_CoordinateAxisType','Time');
% for v = 1:length(Vars_OA)
%     nccreate(filename,Vars_OA{v},'Dimensions',{'Lon',US_RFR_LME.dim.x,'Lat',US_RFR_LME.dim.y,'Time',US_RFR_LME.dim.z});
%     ncwrite(filename,Vars_OA{v},US_RFR_LME.(Vars_OA{v}));
% end
% 
% %% save OA indicator uncertainties as NetCDFs
% filename = ['Data/US_RFR_LME_Inds_Uncer_' date '.nc'];
% time = datenum([repmat(1998,US_RFR_LME.dim.z,1) US_RFR_LME.month+0.5, ...
%     repmat(15,US_RFR_LME.dim.z,1)]);
% nccreate(filename,'Lon','Dimensions',{'Lon',US_RFR_LME.dim.x});
% ncwrite(filename,'Lon',US_RFR_LME.lon);
% ncwriteatt(filename,'Lon','_CoordinateAxisType','Lon');
% nccreate(filename,'Lat','Dimensions',{'Lat',US_RFR_LME.dim.y});
% ncwrite(filename,'Lat',US_RFR_LME.lat);
% ncwriteatt(filename,'Lat','_CoordinateAxisType','Lat');
% nccreate(filename,'Time','Dimensions',{'Time',US_RFR_LME.dim.z});
% ncwrite(filename,'Time',time);
% ncwriteatt(filename,'Time','_CoordinateAxisType','Time');
% for v = 1:length(uVars_OA)
%     nccreate(filename,uVars_OA{v},'Dimensions',{'Lon',US_RFR_LME.dim.x,'Lat',US_RFR_LME.dim.y,'Time',US_RFR_LME.dim.z});
%     ncwrite(filename,uVars_OA{v},US_RFR_LME.(uVars_OA{v}));
% end

%% plot for sanity (requires mapping toolbox)
% figure; hold on;
% worldmap([US_RFR_LME.lim.latmin US_RFR_LME.lim.latmax],...
%     [US_RFR_LME.lim.lonmin US_RFR_LME.lim.lonmax]);
% pcolorm(double(US_RFR_LME.lat),double(US_RFR_LME.lon),...
%     double(mean(US_RFR_LME.fCO2,3,'omitnan'))');
% plot_land('map');
% close

%% Clean up
clear
