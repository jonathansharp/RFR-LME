% Obtain mixed layer depth from GLORYS reanalysis
function Salinity = import_SSS_GLORYS(lat,lon,month,ocean_mask,path)

load([path '/data_to_use/SSS_GLORYS.mat'],'SSS');

% Convert longitude
SSS.lon = convert_lon(SSS.lon);

% Reorder according to longitude
[SSS.lon,lonidx] = sort(SSS.lon);
SSS.sss = SSS.sss(lonidx,:,:);

% Format latitude and longitude
SSS.latitude = repmat(double(SSS.lat)',length(SSS.lon),1);
SSS.longitude = repmat(double(SSS.lon),1,length(SSS.lat));

% Cut down dataset to limits of LME
lonidx = SSS.lon >= min(lon) & SSS.lon <= max(lon);
latidx = SSS.lat >= min(lat) & SSS.lat <= max(lat);
SSS.sss = SSS.sss(lonidx,latidx,:);
SSS.latitude = SSS.latitude(lonidx,latidx,:);
SSS.longitude = SSS.longitude(lonidx,latidx,:);

% Pre-allocate
Salinity = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid
for t = 1:length(month)
    % use climatology for 2021
    if t <= 276
        t_tmp=t;
    else
        t_tmp=t-276:12:t-12;
    end
    % Index where SSS is true
    idx = ~isnan(mean(SSS.sss(:,:,t_tmp),3));
    % Get teporary SSS
    sss_tmp = mean(SSS.sss(:,:,t_tmp),3);
    % Create interpolant over than range
    interp = scatteredInterpolant(SSS.latitude(idx),SSS.longitude(idx),sss_tmp(idx));
    % Fill grid with interpolated SSS
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    sss_tmp = interp(lat_tmp,lon_tmp);
    % Remove values outside of ocean mask
    sss_tmp(~ocean_mask) = NaN;
    Salinity(:,:,t) = sss_tmp;
end

end