% Import ERA5 Winds
function Wind_Speed = import_Winds_ERA5(lat,lon,month,ocean_mask,path)

load([path '/data_to_use/ERA5.mat'],'ERA5');

% Format latitude and longitude
ERA5.latitude = repmat(ERA5.lat',size(ERA5.speed,1),1,size(ERA5.speed,3));
ERA5.longitude = repmat(ERA5.lon,1,size(ERA5.speed,2),size(ERA5.speed,3));

% Cut down dataset to limits of LME
lonidx = ERA5.lon >= min(lon) & ERA5.lon <= max(lon);
latidx = ERA5.lat >= min(lat) & ERA5.lat <= max(lat);
ERA5.speed = ERA5.speed(lonidx,latidx,:);
ERA5.latitude = ERA5.latitude(lonidx,latidx,:);
ERA5.longitude = ERA5.longitude(lonidx,latidx,:);

% Pre-allocate
Wind_Speed = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid
for t = 1:length(ERA5.date)
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),...
        fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.speed(:,:,t)));
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    Wind_Speed(:,:,t) = interp(lon_tmp,lat_tmp);
end

% Remove values outside of ocean mask
for t = 1:length(month)
    speed_tmp = Wind_Speed(:,:,t);
    speed_tmp(~ocean_mask) = NaN;
    Wind_Speed(:,:,t) = speed_tmp;
end

end