% Obtain mixed layer depth from HYCOM model
function Mixed_Layer = import_MLD_HYCOM(lat,lon,month,ocean_mask,path)

load([path '/Data/global/MLD.mat'],'MLD');

% Convert longitude
MLD.lon = convert_lon(MLD.lon);

% Reorder according to longitude
[MLD.lon,lonidx] = sort(MLD.lon);
MLD.mld = MLD.mld(:,lonidx,:);

% Format latitude and longitude
MLD.latitude = repmat(MLD.lat',1,length(MLD.lon));
MLD.longitude = repmat(MLD.lon,length(MLD.lat),1);

% Cut down dataset to limits of LME
lonidx = MLD.lon >= min(lon) & MLD.lon <= max(lon);
latidx = MLD.lat >= min(lat) & MLD.lat <= max(lat);
MLD.mld = MLD.mld(latidx,lonidx,:);
MLD.latitude = MLD.latitude(latidx,lonidx,:);
MLD.longitude = MLD.longitude(latidx,lonidx,:);

% Pre-allocate
Mixed_Layer = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid
for t = 1:size(MLD.mld,3)
    % Index where MLD is true
    idx = ~isnan(MLD.mld(:,:,t));
    % Get teporary MLD
    mld_tmp = MLD.mld(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(MLD.longitude(idx),MLD.latitude(idx),mld_tmp(idx));
    % Fill grid with interpolated SSS
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    mld_tmp = interp(lon_tmp,lat_tmp);
    % Remove values outside of ocean mask
    mld_tmp(~ocean_mask) = NaN;
    Mixed_Layer(:,:,t) = mld_tmp;
end

end