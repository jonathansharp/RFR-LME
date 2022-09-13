% Import ECCO2 SSS
function SSS = import_SSS_ECCO(lat,lon,ocean_mask,path)

load([path '/Data/ECCO_SSS.mat'],'ECCO_SSS');

% Format latitude and longitude
ECCO_SSS.latitude = repmat(ECCO_SSS.lat',size(ECCO_SSS.sss_mon,1),1,size(ECCO_SSS.sss_mon,3));
ECCO_SSS.longitude = repmat(ECCO_SSS.lon,1,size(ECCO_SSS.sss_mon,2),size(ECCO_SSS.sss_mon,3));

% Cut down dataset to limits of LME
lonidx = ECCO_SSS.lon >= min(lon) & ECCO_SSS.lon <= max(lon);
latidx = ECCO_SSS.lat >= min(lat) & ECCO_SSS.lat <= max(lat);
ECCO_SSS.sss_mon = double(ECCO_SSS.sss_mon(lonidx,latidx,:));
ECCO_SSS.latitude = ECCO_SSS.latitude(lonidx,latidx,:);
ECCO_SSS.longitude = ECCO_SSS.longitude(lonidx,latidx,:);

% Interpolate over some gaps in SSS dataset
for t = 1:length(ECCO_SSS.year_mon)
    % Index where ocean mask and SSS are both true
    idx = ~isnan(ECCO_SSS.sss_mon(:,:,t)) & ocean_mask;
    % Get teporary lat, lon, SSS
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    sss_tmp = ECCO_SSS.sss_mon(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),sss_tmp(idx));
    % Index where SSS is nan
    idx = isnan(ECCO_SSS.sss_mon(:,:,t));
    % Fill that area with interpolated SSS
    sss_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    % Remove values outside of ocean mask
    sss_tmp(~ocean_mask) = NaN;
    ECCO_SSS.sss_mon(:,:,t) = sss_tmp;
end

SSS = ECCO_SSS.sss_mon;

end