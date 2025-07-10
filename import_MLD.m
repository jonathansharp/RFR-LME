% import MLD
function data_interp = import_MLD(dpath,vrs,type,lat,lon,time,yr_end,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/MLD_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'CMEMS')
    data_interp = import_MLD_CMEMS(dpath,lat,lon,time,yr_end);
else
    error('Input variable "type" must be "CMEMS"');
end

% save data file
ncsave_3d(['Data/MLD_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'MLD' data_interp 'mixed layer thickness' 'meters'});

else

data_interp = ncread(['Data/MLD_' type '_' vrs '.nc'],'MLD');

end

% create sst animation
if plot_option == 1
    create_animation('MLD',type,time,lat,lon,data_interp,cmocean('tempo'),[0 200],'Mixed Layer Depth','');
    % create_animation('MLD_anom',type,time,lat,lon,data_interp-mean(data_interp,3,'omitnan'),cmocean('balance'),[-2 2],'Mixed Layer Depth Anomaly','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS MLD
function data_interp = import_MLD_CMEMS(dpath,lat,lon,time,yr_end)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_myint_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > conda deactivate

    % file name
    fname = ['cmems_mod_glo_phy_my_0.083deg_P1M-m_multi-vars_180.00W-' ...
        '179.92E_80.00S-90.00N_0.49m_1998-01-01-2021-06-01.nc'];
    fname_int = ['cmems_mod_glo_phy_myint_0.083deg_P1M-m_multi-vars_180.00W-' ...
        '179.92E_80.00S-90.00N_0.49m_2021-07-01-2024-12-01.nc'];
    
    % load dimensions
    inf = ncinfo([dpath 'CMEMS/' fname]);
    data_lat = ncread([dpath 'CMEMS/' fname],'latitude'); % degrees north
    data_lon = ncread([dpath 'CMEMS/' fname],'longitude'); % degrees east
    data_lon_neg = data_lon; data_lon_pos = data_lon;
    data_lon_neg(data_lon_neg >= 0) = NaN;
    data_lon_pos(data_lon_pos < 0) = NaN;
    data_lon_neg = convert_lon(data_lon_neg);
    data_time = ncread([dpath 'CMEMS/' fname],'time'); % hours since 1950-01-01
    data_time = datenum(1950,1,1,data_time,0,0) + 14; % add 14 days for mid-month
    data_time_int = ncread([dpath 'CMEMS/' fname_int],'time'); % hours since 1950-01-01
    data_time_int = datenum(1950,1,1,data_time_int,0,0) + 14; % add 14 days for mid-month
    
    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-floor(min(lat))));
    [~,idx_maxlat] = min(abs(data_lat-ceil(max(lat))));
    [~,idx_minlon1] = min(abs(data_lon_neg-floor(min(lon))));
    [~,idx_maxlon1] = min(abs(data_lon_neg-ceil(max(lon))));
    [~,idx_minlon2] = min(abs(data_lon_pos-floor(min(lon))));
    [~,idx_maxlon2] = min(abs(data_lon_pos-ceil(max(lon))));
    [~,idx_mintime] = min(abs(data_time-floor(min(time))));
    [~,idx_maxtime] = min(abs(data_time-ceil(max(time))));
    [~,idx_mintime_int] = min(abs(data_time_int-floor(min(time))));
    [~,idx_maxtime_int] = min(abs(data_time_int-datenum(yr_end,12,15))); % last month of last year
    
    % cut down dimensions
    data_lat = data_lat(idx_minlat:idx_maxlat);
    data_lon = [data_lon_pos(idx_minlon2:idx_maxlon2);data_lon_neg(idx_minlon1:idx_maxlon1)];
    data_time = data_time(idx_mintime:idx_maxtime);
    data_time_int = data_time_int(idx_mintime_int:idx_maxtime_int);
    
    % read in science-quality data in two longitude chunks
    data_neglon = ncread([dpath 'CMEMS/' fname],'mlotst',[idx_minlon1 idx_minlat idx_mintime],...
        [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1+idx_maxtime-idx_mintime]);
    data_poslon = ncread([dpath 'CMEMS/' fname],'mlotst',[idx_minlon2 idx_minlat idx_mintime],...
        [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1+idx_maxtime-idx_mintime]);
    
    % read in interim data in two longitude chunks
    data_neglon_int = ncread([dpath 'CMEMS/' fname_int],'mlotst',[idx_minlon1 idx_minlat idx_mintime_int],...
        [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1+idx_maxtime_int-idx_mintime_int]);
    data_poslon_int = ncread([dpath 'CMEMS/' fname_int],'mlotst',[idx_minlon2 idx_minlat idx_mintime_int],...
        [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1+idx_maxtime_int-idx_mintime_int]);
    
    % combine data into final grid
    data = cat(3,cat(1,data_poslon,data_neglon),cat(1,data_poslon_int,data_neglon_int));
    clear data_poslon data_neglon data_poslon_int data_neglon_int
    
    % interpolate onto quarter degree grid
    data_interp = nan(length(lon),length(lat),(yr_end-1997)*12);
    [data_lon_grid,data_lat_grid] = ndgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = ndgrid(lon,lat);
    for t = 1:(yr_end-1997)*12
        data_interp(:,:,t) = griddata(double(data_lon_grid),...
            double(data_lat_grid),double(data(:,:,t)),lon_grid,lat_grid);
    end

end

end
