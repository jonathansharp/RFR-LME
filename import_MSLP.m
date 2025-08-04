% import MSLP
function data_interp = import_MSLP(dpath,vrs,type,lat,lon,time,yr_end,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/MSLP_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'NCEP')
    data_interp = import_MSLP_NCEP(dpath,lat,lon,time,yr_end);
else
    error('Input variable "type" must be "NCEP"');
end

% save data file
ncsave_3d(['Data/MSLP_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'MSLP' data_interp 'mean sea level pressure' 'atmospheres'});

else

data_interp = ncread(['Data/MSLP_' type '_' vrs '.nc'],'MSLP');

end

% create sst animation
if plot_option == 1
    create_animation('MSLP',type,time,lat,lon,data_interp,cmocean('rain'),[.9 1.1],'Mean Sea Level Pressure','');
    % create_animation('MSLP_anom',type,time,lat,lon,data_interp-mean(data_interp,3,'omitnan'),cmocean('balance'),[-1 1],'Mean Sea Level Pressure Anomaly','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import NCEP MSLP
function data_interp = import_MSLP_NCEP(dpath,lat,lon,time,yr_end)
    % files obtained from:
    % https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/surface/

    % file name
    fname = ['mslp.mon.mean.nc'];
    
    % load dimensions
    inf = ncinfo([dpath 'NCEP-DOE/' fname]);
    data_lat = ncread([dpath 'NCEP-DOE/' fname],'lat'); % degrees north
    data_lon = ncread([dpath 'NCEP-DOE/' fname],'lon'); % degrees east
    data_time = ncread([dpath 'NCEP-DOE/' fname],'time'); % hours since 1950-01-01
    data_time = datenum(1800,1,1,data_time,0,0) + 14; % add 14 days for mid-month
    
    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-floor(min(lat))));
    [~,idx_maxlat] = min(abs(data_lat-ceil(max(lat))));
    [~,idx_minlon] = min(abs(data_lon-floor(min(lon))));
    [~,idx_maxlon] = min(abs(data_lon-ceil(max(lon))));
    [~,idx_mintime] = min(abs(data_time-floor(min(time))));
    [~,idx_maxtime] = min(abs(data_time-datenum(yr_end,12,15))); % last month of last year

    % cut down dimensions
    data_lat = data_lat(idx_maxlat:idx_minlat);
    data_lon = data_lon(idx_minlon:idx_maxlon);
    % data_time = data_time(idx_mintime:idx_maxtime);
    
    % read in data
    data = ncread([dpath 'NCEP-DOE/' fname],'mslp',[idx_minlon idx_maxlat idx_mintime],...
        [1+idx_maxlon-idx_minlon 1+idx_minlat-idx_maxlat 1+idx_maxtime-idx_mintime]);

    % interpolate onto quarter degree grid
    data_interp = nan(length(lon),length(lat),(yr_end-1997)*12);
    [data_lon_grid,data_lat_grid] = ndgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = ndgrid(lon,lat);
    for t = 1:(yr_end-1997)*12
        data_interp(:,:,t) = griddata(double(data_lon_grid),...
            double(data_lat_grid),double(data(:,:,t)),lon_grid,lat_grid);
    end
    
    % convert Pascals to Atmospheres
    data_interp = data_interp./101325;

end

end
