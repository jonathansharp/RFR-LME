% import SSH
function import_CHL(dpath,vrs,type,lat,lon,time,yr_end,ocean_mask,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/CHL_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'NASA')
    options = weboptions('Username','***','Password','***');
    data_interp = import_CHL_NASA(dpath,lat,lon,time,yr_end,options,ocean_mask);
elseif strcmp(type,'CMEMS')
    data_interp = import_CHL_CMEMS(dpath,lat,lon,time,yr_end,ocean_mask);
else
    error('Input variable "type" must be "CMEMS" or "ECCO"');
end

% save data file
ncsave_3d(['Data/CHL_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'CHL' data_interp 'sea surface chlorophyll' 'milligrams per meter squared'});

else

data_interp = ncread(['Data/CHL_' type '_' vrs '.nc'],'CHL');

end

% create chl animation
if plot_option == 1
    create_animation('CHL',type,time,lat,lon,data_interp,cmocean('algae'),[0 3],'Sea Surface Chlorophyll','mg/m2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import NASA CHL
    function data_interp = import_CHL_NASA(dpath,lat,lon,time,yr_end,options,ocean_mask)

    % Download netcdf files (SeaWiFS)
    fpath = 'Sat/SeaWiFS/Mapped/Monthly/9km/chlor_a/';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3 Mapped/';
    data_lon = ncread([dpath fpath 'S19980011998031.L3m_MO_CHL_chlor_a_9km.nc'],'lon');
    data_lat = ncread([dpath fpath 'S19980011998031.L3m_MO_CHL_chlor_a_9km.nc'],'lat');
    
    % separate longitudes into positive and negative 
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

    % create temporary netcdf
    if isfile('Data/CHL_temp.nc'); delete('Data/CHL_temp.nc'); end
    nccreate('Data/CHL_temp.nc','lon','Dimensions',{'lon',length(data_lon)});
    ncwrite('Data/CHL_temp.nc','lon',data_lon);
    nccreate('Data/CHL_temp.nc','lat','Dimensions',{'lat',length(data_lat)});
    ncwrite('Data/CHL_temp.nc','lat',data_lat);
    nccreate('Data/CHL_temp.nc','time','Dimensions',{'time',Inf});
    nccreate('Data/CHL_temp.nc','chl','Dimensions',{'lon',length(data_lon),...
        'lat',length(data_lat),'time',Inf});
    
    % Define normal and leap year days of year to start month
    year = 1998:2002;
    monthdaynorm = [1 32 60 91 121 152 182 213 244 274 305 335 366];
    monthdayleap = [1 32 61 92 122 153 183 214 245 275 306 336 367];
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
                % this does not work!
                websave([dpath fpath fname],[url fname],options);
            end
            data_temp_neglon = ncread([dpath fpath fname],'chlor_a',...
                [idx_minlon1 idx_maxlat],...
                [1+idx_maxlon1-idx_minlon1 1+idx_minlat-idx_maxlat]);
            data_temp_poslon = ncread([dpath fpath fname],'chlor_a',...
                [idx_minlon2 idx_maxlat],...
                [1+idx_maxlon2-idx_minlon2 1+idx_minlat-idx_maxlat]);
            % append to temporary netcdf
            ncwrite('Data/CHL_temp.nc','chl',[data_temp_poslon;...
                data_temp_neglon],[1 1 (y-1)*12+m]);
        end
    end
    
    % Define month and year limits (MODIS)
    year = 2003:yr_end;
    daynorm = [31 28 31 30 31 30 31 31 30 31 30 31];
    dayleap = [31 29 31 30 31 30 31 31 30 31 30 31];
    % Download netcdf files (MODIS)
    fpath = 'Sat/MODIS/Aqua/Mapped/Monthly/9km/chlor_a/';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://oceandata.sci.gsfc.nasa.gov/getfile/';
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
                % this does not work!
                websave([dpath fpath fname],[url fname],options);
            end
            data_temp_neglon = ncread([dpath fpath fname],'chlor_a',...
                [idx_minlon1 idx_maxlat],...
                [1+idx_maxlon1-idx_minlon1 1+idx_minlat-idx_maxlat]);
            data_temp_poslon = ncread([dpath fpath fname],'chlor_a',...
                [idx_minlon2 idx_maxlat],...
                [1+idx_maxlon2-idx_minlon2 1+idx_minlat-idx_maxlat]);
            % append to temporary netcdf
            ncwrite('Data/CHL_temp.nc','chl',[data_temp_poslon;...
                data_temp_neglon],[1 1 60+(y-1)*12+m]);
        end
    end

    % Pre-allocate
    chl_ts = squeeze(ncread('Data/CHL_temp.nc','chl'));   

    % Interpolate over some gaps in CHL dataset (linear, 1-D, time), then
    % remaining gaps at either end (nearest, 1-D, time)
    for g = 1:size(data_lon)
        for h = 1:size(data_lat)
            if sum(~isnan(chl_ts(g,h,:))) >= 200 % check for "too many" NaNs
                % linear interpolation
                idx = ~isnan(chl_ts(g,h,:));
                chl_tmp = interp1(time(idx),squeeze(chl_ts(g,h,idx)),time,'linear');
                % then, nearest neighbor interpolation
                idx = ~isnan(chl_tmp);
                Chlfit = interp1(time(idx),chl_tmp(idx),time,'nearest','extrap');
                chl_ts(g,h,:) = Chlfit;
            else
                chl_ts(g,h,:) = NaN;
            end
        end
    end

    % adjust ocean mask
    for g = 1:size(lon)
        idx_lon = data_lon >= lon(g)-0.125 & data_lon <= lon(g)+0.125;
        for h = 1:size(lat)
            idx_lat = data_lat >= lat(h)-0.125 & data_lat <= lat(h)+0.125;
            chl_test = chl_ts(idx_lon,idx_lat,1);
            if any(~isnan((chl_test(:))))
                % do nothing
            else
                ocean_mask(g,h) = false;
            end
        end
    end
                
    % Pre-allocate
    data_interp = nan(length(lon),length(lat),(yr_end-1997)*12);
    
    % Interpolate onto quarter degree grid
    for t = 1:(yr_end-1997)*12
        chl_temp = chl_ts(:,:,t);
        % Index where CHL is both true
        idx = ~isnan(chl_temp);
        % Get teporary lat, lon, CHL
        [data_lon2d,data_lat2d] = ndgrid(data_lon,data_lat);
        % Create interpolant over than range
        interp = scatteredInterpolant(data_lon2d(idx),data_lat2d(idx),chl_temp(idx));
        % Get teporary lat and lon
        [lon2d,lat2d] = ndgrid(lon,lat);
        % interpolate to grid
        chl_temp = interp(lon2d,lat2d);
        % remove values outside of ocean mask
        chl_temp(~ocean_mask) = NaN;
        data_interp(:,:,t) = chl_temp;
    end

    % delete temporary file
    delete('Data/CHL_temp.nc');
        
end

% embedded function to import CMEMS CHL
function data_interp = import_CHL_CMEMS(dpath,lat,lon,time,yr_end,ocean_mask)
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
    data_interp = nan(length(lon),length(lat),(yr_end-1997)*12);
    [data_lon_grid,data_lat_grid] = ndgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = ndgrid(lon,lat);

    % read in data
    date = datevec(time);
    year = date(:,1);
    month = date(:,2);
    for t = 1:(yr_end-1997)*12

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
                [idx_maxlon1-idx_minlon1 idx_minlat-idx_maxlat 1]);
            data_poslon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon2 idx_maxlat 1],...
                [idx_maxlon2-idx_minlon2 idx_minlat-idx_maxlat 1]);
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
                [idx_maxlon1-idx_minlon1 idx_minlat-idx_maxlat 1]);
            data_poslon = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_myint_l4-multi-4km_P1M.nc'], ...
                'CHL',[idx_minlon2 idx_maxlat 1],...
                [idx_maxlon2-idx_minlon2 idx_minlat-idx_maxlat 1]);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/' num2str(year(t)) sprintf('%02d',month(t)) '01-' ...
                num2str(year(t)) sprintf('%02d',month(t)) ...
                num2str(day(t-12*(year(t)-1998))) '_' ...
                'cmems_obs-oc_glo_bgc-plankton_myint_l4-multi-4km_P1M.nc'], ...
                'time');
        end
        % interpolate onto quarter degree grid
%         data_interp(:,:,t) = griddata(double(data_lon_grid),...
%             double(data_lat_grid),double([data_poslon;data_neglon]),...
%             lon_grid,lat_grid);
        data_temp = [data_poslon;data_neglon];
        idx = ~isnan(data_temp);
        interp_temp = scatteredInterpolant(double(data_lon_grid(idx)),...
            double(data_lat_grid(idx)),double(data_temp(idx)));
        data_interp_temp = interp_temp(lon_grid,lat_grid);
        data_interp_temp(~ocean_mask) = NaN;
    end

end

end
