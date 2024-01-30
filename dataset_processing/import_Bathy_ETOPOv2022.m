% Obtain bathymetry from ETOPO2
function Bathy = import_Bathy_ETOPOv2022(lat,lon,ocean_mask,path)

ETOPOv2022.lon = ncread([path '/data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc'],'lon');
ETOPOv2022.lat = ncread([path '/data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc'],'lat');
ETOPOv2022.bottomdepth = ncread([path '/data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc'],'z');

% Convert longitude
ETOPOv2022.lon = convert_lon(ETOPOv2022.lon);

% Reorder according to longitude
[ETOPOv2022.lon,lonidx] = sort(ETOPOv2022.lon);
ETOPOv2022.bottomdepth = ETOPOv2022.bottomdepth(lonidx,:);

% Format latitude and longitude
ETOPOv2022.latitude = repmat(ETOPOv2022.lat',length(ETOPOv2022.lon),1);
ETOPOv2022.longitude = repmat(ETOPOv2022.lon,1,length(ETOPOv2022.lat));

% Cut down dataset to limits of LME
lonidx = ETOPOv2022.lon >= min(lon) & ETOPOv2022.lon <= max(lon);
latidx = ETOPOv2022.lat >= min(lat) & ETOPOv2022.lat <= max(lat);
ETOPOv2022.bottomdepth = ETOPOv2022.bottomdepth(lonidx,latidx,:);
ETOPOv2022.latitude = ETOPOv2022.latitude(lonidx,latidx,:);
ETOPOv2022.longitude = ETOPOv2022.longitude(lonidx,latidx,:);

% Interpolate onto quarter degree grid
interp = griddedInterpolant(ETOPOv2022.longitude,...
    ETOPOv2022.latitude,ETOPOv2022.bottomdepth);
lon_tmp = repmat(lon,1,length(lat));
lat_tmp = repmat(lat',length(lon),1);
Bathy = -interp(lon_tmp,lat_tmp);

% Remove values outside of ocean mask
Bathy(~ocean_mask) = NaN;

end
