% import Bathymetry
function data_interp = import_Bathy(dpath,vrs,type,lat,lon,time,varargin)

if ~isfile(['Data/Bathy_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'ETOPO')
    data_interp = import_Bathy_ETOPO(dpath,lat,lon);
else
    error('Input variable "type" must be "ETOPO"');
end

% save data file
ncsave_2d(['Data/Bathy_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},{'Bathy' data_interp 'seafloor depth' 'meters'});

else

data_interp = ncread(['Data/Bathy_' type '_' vrs '.nc'],'Bathy');

end

% plot bathymetry
create_map('Bathy',type,lat,lon,data_interp,cmocean('deep'),[0 6000],'Bottom Depth','m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import ETOPO Bathymetry
    function data_interp = import_Bathy_ETOPO(dpath,lat,lon)
    % file obtained from:
    % https://www.ngdc.noaa.gov/thredds/catalog/global/ETOPO2022/60s/60s_bed_elev_netcdf/catalog.html?dataset=globalDatasetScan/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc

    data_lon = ncread([dpath '/ETOPO/ETOPO_2022_v1_60s_N90W180_bed.nc'],'lon');
    data_lat = ncread([dpath '/ETOPO/ETOPO_2022_v1_60s_N90W180_bed.nc'],'lat');
    data = ncread([dpath '/ETOPO/ETOPO_2022_v1_60s_N90W180_bed.nc'],'z');
    
    % Convert longitude
    data_lon = convert_lon(data_lon,'0-360');
    
    % Reorder according to longitude
    [data_lon,lonidx] = sort(data_lon);
    data = data(lonidx,:);
    
    % Format latitude and longitude
    [data_lon_2d,data_lat_2d] = ndgrid(data_lon,data_lat);
    [lon_2d,lat_2d] = ndgrid(lon,lat);

    % Interpolate onto quarter degree grid
    interp = griddedInterpolant(data_lon_2d,data_lat_2d,data);
    data_interp = -interp(lon_2d,lat_2d);

end

end