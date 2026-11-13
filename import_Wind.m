% import Winds
function data = import_Wind(dpath,vrs,type,lat,lon,time,cds,yr_end,varargin)

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
    data = import_Wind_ERA5(dpath,lat,lon,time,cds,yr_end);
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

% embedded function to import Winds
function data = import_Wind_ERA5(dpath,lat,lon,time,cds,yr_end)

    % Data manually downloaed from:
    % https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview

    % Automated data download not yet working...
    % % set up authentication headers
    % options = weboptions(...
    %     'HeaderFields', {'Authorization', ['Bearer ' cds.token]; 'Content-Type', 'application/json'}, ...
    %     'MediaType', 'application/json', ...
    %     'Timeout', 60);            
    % dataset = 'reanalysis-era5-single-levels';
    % apiUrl = ['https://cds.climate.copernicus.eu/api/v2/resources/' dataset];

    % % define body request structure matching the web API UI
    % years_to_download = arrayfun(@num2str, 1998:yr_end, 'UniformOutput', 0);
    % requestBody = struct('product_type',{{'monthly_averaged_reanalysis'}},...
    %     'variable',{{'10m_wind_speed'}},...
    %     'year',{years_to_download},...
    %     'month',{{'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'}},...
    %     'time',{{'00:00'}},...
    %     'data_format',{{'netcdf'}},...
    %     'download_format',{{'unarchived'}});

    % % submit job request
    % response = webwrite(apiUrl, jsonencode(requestBody), options);
    % statusUrl = response.location; % Get task status URL

    % % poll for completion
    % disp('Processing query on CDS servers...');
    % isComplete = false;
    % while ~isComplete
    %     pause(5); % Wait 5 seconds between checks
    %     status = webread(statusUrl, options);
    %     if strcmp(status.state, 'completed')
    %         downloadUrl = status.location;
    %         isComplete = true;
    %     elseif strcmp(status.state, 'failed')
    %         error('Data request failed on Copernicus server.');
    %     end
    % end

    % % save file locally using websave
    % disp('Downloading Wind data file...');
    % outputFile = 'ERA5_Wind.nc';
    % websave(outputFile, downloadUrl, options);
    % disp('Download finished successfully.');

    % obtain wind data from downloaded file
    data_lat = ncread([dpath 'ERA5_Wind.nc'],'latitude');
    data_lon = ncread([dpath 'ERA5_Wind.nc'],'longitude');
    data_time = ncread([dpath 'ERA5_Wind.nc'],'valid_time');
    data_time = datenum(1970,1,15,0,0,double(data_time));

    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon] = min(abs(data_lon-min(lon)));
    [~,idx_maxlon] = min(abs(data_lon-max(lon)));
    [~,idx_mintime] = min(abs(data_time-min(time)));
    [~,idx_maxtime] = min(abs(data_time-datenum(yr_end,12,15)));
    
    % read in data
    data = ncread([dpath 'ERA5_Wind.nc'],'si10',[idx_minlon idx_maxlat idx_mintime],...
        [1+idx_maxlon-idx_minlon 1+idx_minlat-idx_maxlat 1+idx_maxtime-idx_mintime]);
    data = fliplr(data);
    data(data<-10^6) = NaN; % define NaNs

end

end
