% Convert MATLAB RFR-LME files to NetCDFs

function matlab_to_netcdf(vrs,region)

% load SOCAT grid
load(['Data/' vrs '_gridded'],'SOCAT_grid');

% copy variables to US LME file
US_RFR_LME.lim = SOCAT_grid.lim;
US_RFR_LME.dim = SOCAT_grid.dim;
US_RFR_LME.lon = single(SOCAT_grid.lon);
US_RFR_LME.lat = single(SOCAT_grid.lat);
US_RFR_LME.month = single(SOCAT_grid.month);
US_RFR_LME.year = single(SOCAT_grid.year);
clear SOCAT_grid

% pre-allocate full grid of LME numbers
US_RFR_LME.LME = single(nan(US_RFR_LME.dim.x,US_RFR_LME.dim.y));

% pre-allocate full grid of predictors
Vars_pred = {'SSS' 'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'Bathy' ... 
    'MLD' 'MSLP' 'apCO2'};
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

% pre-allocate full grid of OA Indicators and Uncertainties
% variable names for matlab structure
Vars_OA = {'fCO2' 'pCO2' 'DIC' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
uVars_OA = {'ufCO2' 'upCO2' 'uDIC' 'uTA' 'upH' 'uOmA' 'uOmC' 'uH' 'uCO3' 'uRF'};
% variable names for NetCDF file
varnames_OA = {'fco2' 'pco2' 'dic' 'ta' 'ph' 'om_a' 'om_c' 'h' 'co3' 'rf'};
varnames_uOA = {'fco2_u' 'pco2_u' 'dic_u' 'ta_u' 'ph_u' 'om_a_u' 'om_c_u' ...
    'h_u' 'co3_u' 'rf_u'};
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
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    load(['Data/RFR-LME/' vrs '_' region{n}],'RFR_LME');

    %% define LME mask
    idx_lon = US_RFR_LME.lon >= min(LME.lon) & US_RFR_LME.lon <= max(LME.lon);
    idx_lat = US_RFR_LME.lat >= min(LME.lat) & US_RFR_LME.lat <= max(LME.lat);

    %% Add regional predictors to full grid
    for v = 1:length(Vars_pred)
        tmp_var = US_RFR_LME.(Vars_pred{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            LME.(Vars_pred{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(Vars_pred{v})(idx) = tmp_var(idx);
    end    

    %% Add OA indicators to full grid
    for v = 1:length(Vars_OA)
        tmp_var = US_RFR_LME.(Vars_OA{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            RFR_LME.(Vars_OA{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(Vars_OA{v})(idx) = tmp_var(idx);
    end

    %% Add OA indicator uncertainties to full grid
    for v = 1:length(uVars_OA)
        tmp_var = US_RFR_LME.(uVars_OA{v});
        tmp_var(idx_lon,idx_lat,:) = ...
            RFR_LME.(uVars_OA{v});
        idx = ~isnan(tmp_var);
        US_RFR_LME.(uVars_OA{v})(idx) = tmp_var(idx);
    end

    %% add LME number to grid
    tmp_idx = false(size(US_RFR_LME.LME));
    tmp_idx(idx_lon,idx_lat) = ...
        any(RFR_LME.idxspc,3);
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
% define file name
filename = ['Data/NetCDFs_' date '/US_RFR_LME.nc'];
% delete old file if it exists
if isfile(filename); delete(filename); end
% longitude
nccreate(filename,'lon','Dimensions',{'lon',US_RFR_LME.dim.x},...
    'DataType','single','FillValue',-9999.0);
ncwrite(filename,'lon',US_RFR_LME.lon);
%ncwriteatt(filename,'lon','bounds','lon_bounds');
ncwriteatt(filename,'lon','units','degrees_east');
ncwriteatt(filename,'lon','axis','X');
ncwriteatt(filename,'lon','long_name','longitude');
ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
% latitude
nccreate(filename,'lat','Dimensions',{'lat',US_RFR_LME.dim.y},...
    'DataType','single','FillValue',-9999.0);
ncwrite(filename,'lat',US_RFR_LME.lat);
%ncwriteatt(filename,'lat','bounds','lat_bounds');
ncwriteatt(filename,'lat','units','degrees_north');
ncwriteatt(filename,'lat','axis','Y');
ncwriteatt(filename,'lat','long_name','latitude');
ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
% time
nccreate(filename,'time','Dimensions',{'time',US_RFR_LME.dim.z},...
    'DataType','single','FillValue',-9999.0);
ncwrite(filename,'time',time);
%ncwriteatt(filename,'time','bounds','time_bounds');
ncwriteatt(filename,'time','units','days since 1-Jan-0000');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');
% predictor variable
for v = 1:length(Vars_pred)
    if ~strcmp(Vars_pred{v},'Bathy') && ~strcmp(Vars_pred{v},'Dist')
        nccreate(filename,varnames_pred{v},'Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
                'time',US_RFR_LME.dim.z},'DataType','single','FillValue',-9999.0,...
                    'DeflateLevel',5,'Shuffle',true,...
                        'ChunkSize',[US_RFR_LME.dim.x,US_RFR_LME.dim.y,1]);
    else
        nccreate(filename,varnames_pred{v},'Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y},...
                'DataType','single','FillValue',-9999.0,...
                    'DeflateLevel',5,'Shuffle',true);
    end
    ncwrite(filename,varnames_pred{v},US_RFR_LME.(Vars_pred{v}));
    ncwriteatt(filename,varnames_pred{v},'units',units_pred{v});
    ncwriteatt(filename,varnames_pred{v},'long_name',names_pred{v});
end
for v = 1:length(Vars_OA)
    % variable
    nccreate(filename,varnames_OA{v},'Dimensions',...
        {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
        'time',US_RFR_LME.dim.z},'DataType','single','FillValue',-9999.0,...
        'DeflateLevel',5,'Shuffle',true,...
        'ChunkSize',[US_RFR_LME.dim.x,US_RFR_LME.dim.y,1]);
    ncwrite(filename,varnames_OA{v},US_RFR_LME.(Vars_OA{v}));
    ncwriteatt(filename,varnames_OA{v},'units',units_OA{v});
    ncwriteatt(filename,varnames_OA{v},'long_name',names_OA{v});
    % variable uncertainty
    nccreate(filename,varnames_uOA{v},'Dimensions',...
        {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y,...
        'time',US_RFR_LME.dim.z},'DataType','single','FillValue',-9999.0,...
        'DeflateLevel',5,'Shuffle',true,...
        'ChunkSize',[US_RFR_LME.dim.x,US_RFR_LME.dim.y,1]);
    ncwrite(filename,varnames_uOA{v},US_RFR_LME.(uVars_OA{v}));
    ncwriteatt(filename,varnames_uOA{v},'units',units_OA{v});
    ncwriteatt(filename,varnames_uOA{v},'long_name',['Uncertainty in ' names_OA{v}]);
    % LME number
    nccreate(filename,'lme_num','Dimensions',...
            {'lon',US_RFR_LME.dim.x,'lat',US_RFR_LME.dim.y},...
            'DataType','single','FillValue',-9999.0);
    ncwrite(filename,'lme_num',US_RFR_LME.LME);
    ncwriteatt(filename,'lme_num','reference_table','lme_ref');
    ncwriteatt(filename,'lme_num','long_name','Large Marine Ecosystem Number');
    ncwriteatt(filename,'lme_num','information','See "lme_ref" for the order of LMEs associated with each LME number.');
    % LME reference
    nccreate(filename,'lme_ref','Dimensions',{'Large Marine Ecosystem Names',11,'b',4},'DataType','char');
    ncwrite(filename,'lme_ref',char(region));
    ncwriteatt(filename,'lme_ref','long_name','Large Marine Ecosystem Reference Table');
    ncwriteatt(filename,'lme_ref','information','This variable provides the order of LME abbreviations associated with each LME number.');
end

end