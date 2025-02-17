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

% create ssh animation
if plot_option == 1
    create_animation('CHL',type,time,lat,lon,log10(data),cmocean('algae'),[0.001 10],'Sea Surface Chlorophyll','mg/m2');
end

% save data file
ncsave_3d(['Data/CHL_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},{'time' time-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'CHL' data 'sea surface chlorophyll' 'milligrams per meter squared'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS SSH
function data = import_CHL_NASA(dpath,lat,lon,time)

    % Define normal and leap year days of year to start month
    year = 1998:2002;
    monthdaynorm = [1 32 60 91 121 152 182 213 244 274 305 335 366];
    monthdayleap = [1 32 61 92 122 153 183 214 245 275 306 336 367];
    
    % Download netcdf files (SeaWiFS)
    fpath = 'Sat/SeaWiFS/Mapped/Monthly/9km/chlor_a';
    if ~isfolder([dpath fpath]); mkdir([dpath fpath]); end
    url = 'https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3 Mapped/';
    for y = 1:length(year)
        fname = 
        if ~isfile([dpath fpath fname])
            websave([dpath fpath fname],[url fname]);
        end
    
        % Download netcdf files (SeaWiFS)
        for y = 1:length(year)
            if  year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
                year(y) == 2012 || year(y) == 2016 || year(y) == 2020
                day = monthdayleap;
            else
                day = monthdaynorm;
            end
            for m = 1:12
                CHL.chl(:,:,(y-1)*12+m) = ...
                    ncread(['Data_Import/chlor_a/S' num2str(year(y)) sprintf('%03d',day(m)) ...
                    num2str(year(y)) sprintf('%03d',day(m+1)-1) ...
                    '.L3m_MO_CHL_chlor_a_9km.nc'],'chlor_a');
            end
        end
    end
    
    % Define month and year limits (MODIS)
    year = 2003:2022;
    daynorm = [31 28 31 30 31 30 31 31 30 31 30 31];
    dayleap = [31 29 31 30 31 30 31 31 30 31 30 31];
    
    % Import netcdf files (MODIS)
    for y = 1:length(year)
        if  year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
            year(y) == 2012 || year(y) == 2016 || year(y) == 2020
            day = dayleap;
        else
            day = daynorm;
        end
        for m = 1:12
            CHL.chl(:,:,(y-1)*12+m) = ...
                ncread(['chlor_a/AQUA_MODIS.' num2str(year(y)) sprintf('%02d',m) ...
                '01_' num2str(year(y)) sprintf('%02d',m) sprintf('%02d',day(m)) ...
                '.L3m.MO.CHL.chlor_a.9km.nc'],'chlor_a');
        end
    end
    % Clean up
    clear m y path dayleap daynorm day year
    
    % Import lat, lon, time
    CHL.lat = ncread(['chlor_a/AQUA_MODIS.20030101_20030131.L3m.MO.CHL.chlor_a.9km.nc'],'lat');
    CHL.lon = ncread(['chlor_a/AQUA_MODIS.20030101_20030131.L3m.MO.CHL.chlor_a.9km.nc'],'lon');
    CHL.time = datenum([repelem(1998:2022,1,12)' repmat(1:12,1,25)' repmat(15,1,25*12)']);

end

% embedded function to import ECCO SSH
function data = import_CHL_CMEMS(dpath,lat,lon,time)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m
    %     > conda deactivate

    % file path
    fpath = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_202112/';
    fpath_int = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m_202311/';

    % load dimensions
    temp_file = '1993/dt_global_allsat_msla_h_y1993_m01.nc';
    inf = ncinfo([dpath fpath temp_file]);
    data_lat = ncread([dpath fpath temp_file],'latitude'); % degrees north
    data_lon = ncread([dpath fpath temp_file],'longitude'); % degrees east
    
    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon] = min(abs(data_lon-min(lon)));
    [~,idx_maxlon] = min(abs(data_lon-max(lon)));

    % cut down dimensions
    data_lat = data_lat(idx_minlat:idx_maxlat);
    data_lon = data_lon(idx_minlon:idx_maxlon);

    % read in data
    date = datevec(time);
    year = date(:,1);
    month = date(:,2);
    for t = 1:length(time)
        try
            % science-quality data
            data(:,:,t) = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                [idx_minlon idx_minlat 1],...
                [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1]);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                '_m' sprintf('%02d',month(t)) '.nc'],'time');
        catch
            try
                % near-real-time data
                data(:,:,t) = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                    [idx_minlon idx_minlat 1],...
                    [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1]);
                data_time(t) = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'time');
            catch
            end
        end
    end

end

end
