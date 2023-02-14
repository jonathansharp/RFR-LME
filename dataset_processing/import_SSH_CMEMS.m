% Import CMEMS SSH
function SSH = import_SSH_CMEMS(lat,lon,ocean_mask,path)

load([path '/Data/global/CMEMS_SSH.mat'],'CMEMS_SSH');

% Format latitude and longitude
CMEMS_SSH.latitude = repmat(CMEMS_SSH.lat',size(CMEMS_SSH.sla_mon,1),1,size(CMEMS_SSH.sla_mon,3));
CMEMS_SSH.longitude = repmat(CMEMS_SSH.lon,1,size(CMEMS_SSH.sla_mon,2),size(CMEMS_SSH.sla_mon,3));

% Cut down dataset to limits of LME
lonidx = CMEMS_SSH.lon >= min(lon) & CMEMS_SSH.lon <= max(lon);
latidx = CMEMS_SSH.lat >= min(lat) & CMEMS_SSH.lat <= max(lat);
CMEMS_SSH.sla_mon = double(CMEMS_SSH.sla_mon(lonidx,latidx,:));
CMEMS_SSH.latitude = CMEMS_SSH.latitude(lonidx,latidx,:);
CMEMS_SSH.longitude = CMEMS_SSH.longitude(lonidx,latidx,:);

% Interpolate over some gaps in SSH dataset
for t = 1:length(CMEMS_SSH.year_mon)
    % Index where ocean mask and SSH are both true
    idx = ~isnan(CMEMS_SSH.sla_mon(:,:,t)) & ocean_mask;
    % Get teporary lat, lon, SSH
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    ssh_tmp = CMEMS_SSH.sla_mon(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),ssh_tmp(idx));
    % Index where SSH is nan
    idx = isnan(CMEMS_SSH.sla_mon(:,:,t));
    % Fill that area with interpolated SSH
    ssh_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    % Remove values outside of ocean mask
    ssh_tmp(~ocean_mask) = NaN;
    CMEMS_SSH.sla_mon(:,:,t) = ssh_tmp;
end

SSH = CMEMS_SSH.sla_mon;

end