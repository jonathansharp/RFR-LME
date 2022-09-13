% Import ECCO2 SSH
function SSH = import_SSH_ECCO(lat,lon,ocean_mask,path)

load([path '/Data/ECCO_SSH.mat'],'ECCO_SSH');

% Format latitude and longitude
ECCO_SSH.latitude = repmat(ECCO_SSH.lat',size(ECCO_SSH.ssh_mon,1),1,size(ECCO_SSH.ssh_mon,3));
ECCO_SSH.longitude = repmat(ECCO_SSH.lon,1,size(ECCO_SSH.ssh_mon,2),size(ECCO_SSH.ssh_mon,3));

% Cut down dataset to limits of LME
lonidx = ECCO_SSH.lon >= min(lon) & ECCO_SSH.lon <= max(lon);
latidx = ECCO_SSH.lat >= min(lat) & ECCO_SSH.lat <= max(lat);
ECCO_SSH.ssh_mon = double(ECCO_SSH.ssh_mon(lonidx,latidx,:));
ECCO_SSH.latitude = ECCO_SSH.latitude(lonidx,latidx,:);
ECCO_SSH.longitude = ECCO_SSH.longitude(lonidx,latidx,:);

% Interpolate over some gaps in SSH dataset
for t = 1:length(ECCO_SSH.year_mon)
    % Index where ocean mask and SSH are both true
    idx = ~isnan(ECCO_SSH.ssh_mon(:,:,t)) & ocean_mask;
    % Get teporary lat, lon, SSH
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    ssh_tmp = ECCO_SSH.ssh_mon(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),ssh_tmp(idx));
    % Index where SSH is nan
    idx = isnan(ECCO_SSH.ssh_mon(:,:,t));
    % Fill that area with interpolated SSH
    ssh_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    % Remove values outside of ocean mask
    ssh_tmp(~ocean_mask) = NaN;
    ECCO_SSH.ssh_mon(:,:,t) = ssh_tmp;
end

SSH = ECCO_SSH.ssh_mon;

end