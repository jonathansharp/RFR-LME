% Import OISSTv2
function SST = import_SST_OISST(lat,lon,ocean_mask)

load('Data/OISSTv2_monthly_averaged_data.mat','OISST');

% Format latitude and longitude
OISST.latitude = repmat(OISST.lat',size(OISST.sst_mon,1),1,size(OISST.sst_mon,3));
OISST.longitude = repmat(OISST.lon,1,size(OISST.sst_mon,2),size(OISST.sst_mon,3));

% Cut down dataset to limits of LME
lonidx = OISST.lon >= min(lon) & OISST.lon <= max(lon);
latidx = OISST.lat >= min(lat) & OISST.lat <= max(lat);
OISST.sst_mon = double(OISST.sst_mon(lonidx,latidx,:));
OISST.latitude = OISST.latitude(lonidx,latidx,:);
OISST.longitude = OISST.longitude(lonidx,latidx,:);

% Interpolate over some gaps in SST dataset
for t = 1:length(OISST.year_mon)
    % Index where ocean mask and SST are both true
    idx = ~isnan(OISST.sst_mon(:,:,t)) & ocean_mask;
    % Get teporary lat, lon, SST
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    sst_tmp = OISST.sst_mon(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),sst_tmp(idx));
    % Index where SST is nan but ocean mask is true
    idx = isnan(OISST.sst_mon(:,:,t)) & ocean_mask;
    % Fill that area with interpolated SST
    sst_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    OISST.sst_mon(:,:,t) = sst_tmp;
end

SST = OISST.sst_mon;

end