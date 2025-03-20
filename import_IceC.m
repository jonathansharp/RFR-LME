% import IceC
function data = import_IceC(dpath,vrs,type,lat,lon,time,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% Import based on "type"
if strcmp(type,'OISST')
    data = import_IceC_OISST(dpath,lat,lon,time);
else
    error('Input variable "type" must be "OISST"');
end

% save data file as netcdf
ncsave_3d(['Data/IceC_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},{'time' time-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'IceC' data 'ice coverage' 'percent'});

% create ice coverage animation
if plot_option == 1
    create_animation('IceC',type,time,lat,lon,data*100,...
        flipud(cmocean('ice')),[0 100],'Ice Coverage','%');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import OISST ice coverage
function data = import_IceC_OISST(dpath,lat,lon,time)

    % obtain OISST ice coverage file if downloaded file is older than one month
    fpath = 'OISST/IceC/';
    fname = 'icec.mon.mean.nc';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://psl.noaa.gov/thredds/fileServer/Datasets/noaa.oisst.v2.highres/';
    if isfile([dpath fpath fname])
        inf = dir([dpath fpath fname]);
        if datenum(inf.date) - datenum(date) > 30
            websave([dpath fpath fname],[url fname]);
        end
    else
        websave([dpath 'OISST/IceC/' fname],[url fname]);
    end

    % load dimensions
    inf = ncinfo([dpath fpath fname]);
    data_lat = ncread([dpath fpath fname],'lat'); % degrees north
    data_lon = ncread([dpath fpath fname],'lon'); % degrees east
    data_time = ncread([dpath fpath fname],'time'); % days since 1800-01-00
    data_time = datenum(1800,1,data_time,0,0,0) + 15; % add 15 days for mid-month

    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon] = min(abs(data_lon-min(lon)));
    [~,idx_maxlon] = min(abs(data_lon-max(lon)));
    [~,idx_mintime] = min(abs(data_time-min(time)));
    [~,idx_maxtime] = min(abs(data_time-max(time)));
    
    % read in data
    data = ncread([dpath fpath fname],'icec',[idx_minlon idx_minlat idx_mintime],...
        [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1+idx_maxtime-idx_mintime]);
    data(data<-10^6) = NaN; % define NaNs

end

end
