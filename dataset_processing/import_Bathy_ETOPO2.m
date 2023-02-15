% Obtain bathymetry from ETOPO2
function Bathy = import_Bathy_ETOPO2(lat,lon,ocean_mask,path)

load([path '/data_to_use/ETOPO2.mat'],'ETOPO2');

% Convert longitude
ETOPO2.lon = convert_lon(ETOPO2.lon);

% Reorder according to longitude
[ETOPO2.lon,lonidx] = sort(ETOPO2.lon);
ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,:);

% Format latitude and longitude
ETOPO2.latitude = repmat(ETOPO2.lat',length(ETOPO2.lon),1);
ETOPO2.longitude = repmat(ETOPO2.lon,1,length(ETOPO2.lat));

% Cut down dataset to limits of LME
lonidx = ETOPO2.lon >= min(lon) & ETOPO2.lon <= max(lon);
latidx = ETOPO2.lat >= min(lat) & ETOPO2.lat <= max(lat);
ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx,:);
ETOPO2.latitude = ETOPO2.latitude(lonidx,latidx,:);
ETOPO2.longitude = ETOPO2.longitude(lonidx,latidx,:);

% Interpolate onto quarter degree grid
interp = griddedInterpolant(ETOPO2.longitude,ETOPO2.latitude,ETOPO2.bottomdepth);
lon_tmp = repmat(lon,1,length(lat));
lat_tmp = repmat(lat',length(lon),1);
Bathy = interp(lon_tmp,lat_tmp);

% Remove values outside of ocean mask
Bathy(~ocean_mask) = NaN;

end
