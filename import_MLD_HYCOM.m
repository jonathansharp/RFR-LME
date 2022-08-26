% Obtain mixed layer depth from HYCOM model
function Mixed_Layer = import_MLD_HYCOM(lat,lon,month,ocean_mask)

load('Data/MLD.mat','MLD');

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
    interp = griddedInterpolant(flipud(MLD.longitude)',...
        flipud(MLD.latitude)',flipud(MLD.mld(:,:,t))');
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    Mixed_Layer(:,:,t) = interp(lon_tmp,lat_tmp);
end
