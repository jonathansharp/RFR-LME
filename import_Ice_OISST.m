% Import Ice Concentration
function IceC = import_Ice_OISST(lat,lon,ocean_mask,path)

load([path '/Data/global/OISST_ice.mat'],'OISST_ice');

% Format latitude and longitude
OISST_ice.latitude = repmat(OISST_ice.lat',size(OISST_ice.icec_mon,1),1,size(OISST_ice.icec_mon,3));
OISST_ice.longitude = repmat(OISST_ice.lon,1,size(OISST_ice.icec_mon,2),size(OISST_ice.icec_mon,3));

% Cut down dataset to limits of LME
lonidx = OISST_ice.lon >= min(lon) & OISST_ice.lon <= max(lon);
latidx = OISST_ice.lat >= min(lat) & OISST_ice.lat <= max(lat);
OISST_ice.icec_mon = double(OISST_ice.icec_mon(lonidx,latidx,:));
OISST_ice.latitude = OISST_ice.latitude(lonidx,latidx,:);
OISST_ice.longitude = OISST_ice.longitude(lonidx,latidx,:);

% Interpolate over some gaps in Ice Concentration dataset
for t = 1:length(OISST_ice.year_mon)
    % Index where ocean mask and Ice Concentration are both true
    idx = ~isnan(OISST_ice.icec_mon(:,:,t)) & ocean_mask;
    % Get teporary lat, lon, Ice Concentration
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    icec_tmp = OISST_ice.icec_mon(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),icec_tmp(idx));
    % Index where Ice Concentration is nan but ocean mask is true
    idx = isnan(OISST_ice.icec_mon(:,:,t));
    % Fill that area with interpolated Ice Concentration
    icec_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));    
    % Remove values outside of ocean mask
    icec_tmp(~ocean_mask) = NaN;
    OISST_ice.icec_mon(:,:,t) = icec_tmp;
end

IceC = OISST_ice.icec_mon;

end