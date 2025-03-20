% import SSH
function data = import_CHL(dpath,vrs,type,lat,lon,time,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% Import based on "type"
if strcmp(type,'NASA')
    data = import_CHL_NASA(dpath,lat,lon,time);
elseif strcmp(type,'CMEMS')
    data = import_CHL_CMEMS(dpath,lat,lon,time);
else
    error('Input variable "type" must be "CMEMS" or "ECCO"');
end

% save data file
ncsave_3d(['Data/CHL_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},{'time' time-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'CHL' data 'sea surface chlorophyll' 'milligrams per meter squared'});

% create chl animation
if plot_option == 1
    create_animation('CHL',type,time,lat,lon,log10(data),cmocean('algae'),[0.001 10],'Sea Surface Chlorophyll','mg/m2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import NASA CHL
function data = import_CHL_NASA(dpath,lat,lon,time)

    % Define normal and leap year days of year to start month
    year = 1998:2002;
    monthdaynorm = [1 32 60 91 121 152 182 213 244 274 305 335 366];
    monthdayleap = [1 32 61 92 122 153 183 214 245 275 306 336 367];
    % Download netcdf files (SeaWiFS)
    fpath = 'Sat/SeaWiFS/Mapped/Monthly/9km/chlor_a/';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3 Mapped/';
    for y = 1:length(year)
        if  year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
            year(y) == 2012 || year(y) == 2016 || year(y) == 2020
            day = monthdayleap;
        else
            day = monthdaynorm;
        end
        for m = 1:12
            fname = ['S' num2str(year(y)) sprintf('%03d',day(m)) ...
                num2str(year(y)) sprintf('%03d',day(m+1)-1) ...
                '.L3m_MO_CHL_chlor_a_9km.nc'];
            if ~isfile([dpath fpath fname])
                websave([dpath fpath fname],[url fname]);
            end
            data(:,:,(y-1)*12+m) = ...
                ncread([dpath fpath fname],'chlor_a');
        end
    end
    
    % Define month and year limits (MODIS)
    year = 2003:2024;
    daynorm = [31 28 31 30 31 30 31 31 30 31 30 31];
    dayleap = [31 29 31 30 31 30 31 31 30 31 30 31];
    % Download netcdf files (MODIS)
    fpath = 'Sat/MODIS/Aqua/Mapped/Monthly/9km/chlor_a/';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3 Mapped/';
    for y = 1:length(year)
        if  year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
            year(y) == 2012 || year(y) == 2016 || year(y) == 2020 || year(y) == 2024
            day = dayleap;
        else
            day = daynorm;
        end
        for m = 1:12
            fname = ['AQUA_MODIS.' num2str(year(y)) sprintf('%02d',m) '01_' ...
                num2str(year(y)) sprintf('%02d',m) sprintf('%02d',day(m)) ...
                '.L3m.MO.CHL.chlor_a.9km.nc'];
            if ~isfile([dpath fpath fname])
                websave([dpath fpath fname],[url fname]);
            end
            data(:,:,60+(y-1)*12+m) = ...
                ncread([dpath fpath fname],'chlor_a');
        end
    end
    % Clean up
    clear m y path dayleap daynorm day year
    
    % Import lat, lon, time
    CHL.lat = ncread(['chlor_a/AQUA_MODIS.20030101_20030131.L3m.MO.CHL.chlor_a.9km.nc'],'lat');
    CHL.lon = ncread(['chlor_a/AQUA_MODIS.20030101_20030131.L3m.MO.CHL.chlor_a.9km.nc'],'lon');
    CHL.time = datenum([repelem(1998:2022,1,12)' repmat(1:12,1,25)' repmat(15,1,25*12)']);



%     % load dimensions
%     temp_file = '1998/19980101_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc';
%     inf = ncinfo([dpath fpath temp_file]);
%     data_lat = ncread([dpath fpath temp_file],'latitude'); % degrees north
%     data_lon = ncread([dpath fpath temp_file],'longitude'); % degrees east
%     
%     % index based on dimensions
%     [~,idx_minlat] = min(abs(data_lat-min(lat)));
%     [~,idx_maxlat] = min(abs(data_lat-max(lat)));
%     [~,idx_minlon] = min(abs(data_lon-min(lon)));
%     [~,idx_maxlon] = min(abs(data_lon-max(lon)));
% 
%     % cut down dimensions
%     data_lat = data_lat(idx_minlat:idx_maxlat);
%     data_lon = data_lon(idx_minlon:idx_maxlon);
% 
%     % read in data
%     date = datevec(time);
%     year = date(:,1);
%     month = date(:,2);

end

% embedded function to import CMEMS CHL
function data = import_CHL_CMEMS(dpath,lat,lon,time)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine get --dataset-id cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M
    %     > copernicusmarine get --dataset-id cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M
    %     > conda deactivate

    % file path
    fpath = 'CMEMS/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_202411/';

    % load dimensions
    temp_file = '1998/19980101-19980131_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc';
    inf = ncinfo([dpath fpath temp_file]);
    data_lat = ncread([dpath fpath temp_file],'lat'); % degrees north
    data_lon = ncread([dpath fpath temp_file],'lon'); % degrees east
    data_lon_neg = data_lon; data_lon_pos = data_lon;
    data_lon_neg(data_lon_neg >= 0) = NaN;
    data_lon_pos(data_lon_pos < 0) = NaN;
    data_lon_neg = convert_lon(data_lon_neg);

    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon1] = min(abs(data_lon_neg-floor(min(lon))));
    [~,idx_maxlon1] = min(abs(data_lon_neg-ceil(max(lon))));
    [~,idx_minlon2] = min(abs(data_lon_pos-floor(min(lon))));
    [~,idx_maxlon2] = min(abs(data_lon_pos-ceil(max(lon))));

    % cut down dimensions
    data_lat = data_lat(idx_maxlat:idx_minlat);
    data_lon = [data_lon_pos(idx_minlon2:idx_maxlon2);data_lon_neg(idx_minlon1:idx_maxlon1)];

    % days of months
    daynorm = [31 28 31 30 31 30 31 31 30 31 30 31];
    dayleap = [31 29 31 30 31 30 31 31 30 31 30 31];

    % create grid for interpolation
    data_interp = nan(length(lon),length(lat),length(time));
    [data_lon_grid,data_lat_grid] = ndgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = ndgrid(lon,lat);

    % read in data
    date = datevec(time);
    year = date(:,1);
    month = date(:,2);
    for t = 1:length(time)

        if  year(t) == 2000 || year(t) == 2004 || year(t) == 2008 || ...
            year(t) == 2012 || year(t) == 2016 || year(t) == 2020 || ...
            year(t) == 2024
            day = dayleap;
        else
            day = daynorm;
        end
        try
            % science-quality data (in two longitude chunks)
            data_neglon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon1 idx_maxlat 1],...
                [1+idx_maxlon1-idx_minlon1 1+idx_minlat-idx_maxlat 1]);
            data_poslon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon2 idx_maxlat 1],...
                [1+idx_maxlon2-idx_minlon2 1+idx_minlat-idx_maxlat 1]);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc'], ...
                'time');
        catch
            % near-real-time data
            data_neglon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_myint_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon1 idx_maxlat 1],...
                [1+idx_maxlon1-idx_minlon1 1+idx_minlat-idx_maxlat 1]);
            data_poslon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_myint_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon2 idx_maxlat 1],...
                [1+idx_maxlon2-idx_minlon2 1+idx_minlat-idx_maxlat 1]);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_myint_l4-multi-4km_P1M.nc'], ...
                'time');
        end
        % interpolate onto quarter degree grid
        data_interp(:,:,t) = griddata(double(data_lon_grid),...
            double(data_lat_grid),double([data_poslon;data_neglon]),lon_grid,lat_grid);
    end

end

end
