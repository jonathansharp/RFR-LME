% import Winds
function data = import_Wind(dpath,vrs,type,lat,lon,time,yr_end,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/Wind_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'ERA5')
    data = import_Wind_ERA5(dpath,lat,lon,time,yr_end);
else
    error('Input variable "type" must be "ERA5"');
end

% save data file
ncsave_3d(['Data/Wind_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'Wind' data 'wind speed' 'meters per second'});

else

data = ncread(['Data/Wind_' type '_' vrs '.nc'],'Wind');

end

% create sst animation
if plot_option == 1
    create_animation('Wind',type,time,lat,lon,data,cmocean('speed'),[-5 35],'Wind Speed','m/s');
    create_animation('Wind_anom',type,time,lat,lon,data-mean(data,3,'omitnan'),cmocean('balance'),[-5 5],'Wind Speed Anomaly','m/s');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import OISST
function data = import_Wind_ERA5(dpath,lat,lon,time,yr_end)

    % obtain wind data
    data_lat = ncread([dpath '/ERA5_2024/ERA5_Wind.nc'],'latitude');
    data_lon = ncread([dpath '/ERA5_2024/ERA5_Wind.nc'],'longitude');
    data_time = ncread([dpath '/ERA5_2024/ERA5_Wind.nc'],'valid_time');
    data_time = datenum(1970,1,15,0,0,double(data_time));

    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon] = min(abs(data_lon-min(lon)));
    [~,idx_maxlon] = min(abs(data_lon-max(lon)));
    [~,idx_mintime] = min(abs(data_time-min(time)));
    [~,idx_maxtime] = min(abs(data_time-datenum(yr_end,12,15)));
    
    % read in data
    data = ncread([dpath '/ERA5_2024/ERA5_Wind.nc'],'si10',[idx_minlon idx_maxlat idx_mintime],...
        [1+idx_maxlon-idx_minlon 1+idx_minlat-idx_maxlat 1+idx_maxtime-idx_mintime]);
    data = fliplr(data);
    data(data<-10^6) = NaN; % define NaNs

end

end
