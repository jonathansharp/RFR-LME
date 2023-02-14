% Obtain mixed layer depth from GLORYS reanalysis
function Mixed_Layer = import_MLD_GLORYS(lat,lon,month,ocean_mask,path)

load([path '/Data/global/MLD_GLORYS.mat'],'MLD');

% Convert longitude
MLD.lon = convert_lon(MLD.lon);

% Reorder according to longitude
[MLD.lon,lonidx] = sort(MLD.lon);
MLD.mld = MLD.mld(lonidx,:,:);

% Format latitude and longitude
MLD.latitude = repmat(double(MLD.lat)',length(MLD.lon),1);
MLD.longitude = repmat(double(MLD.lon),1,length(MLD.lat));

% Cut down dataset to limits of LME
lonidx = MLD.lon >= min(lon) & MLD.lon <= max(lon);
latidx = MLD.lat >= min(lat) & MLD.lat <= max(lat);
MLD.mld = MLD.mld(lonidx,latidx,:);
MLD.latitude = MLD.latitude(lonidx,latidx,:);
MLD.longitude = MLD.longitude(lonidx,latidx,:);

% Pre-allocate
Mixed_Layer = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid
for t = 1:length(month)
    % use climatology for 2021
    if t <= 276
        t_tmp=t;
    else
        t_tmp=t-276:12:t-12;
    end
    % Index where MLD is true
    idx = ~isnan(mean(MLD.mld(:,:,t_tmp),3));
    % Get teporary MLD
    mld_tmp = mean(MLD.mld(:,:,t_tmp),3);
    % Create interpolant over than range
    interp = scatteredInterpolant(MLD.latitude(idx),MLD.longitude(idx),mld_tmp(idx));
    % Fill grid with interpolated SSS
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    mld_tmp = interp(lat_tmp,lon_tmp);
    % Remove values outside of ocean mask
    mld_tmp(~ocean_mask) = NaN;
    Mixed_Layer(:,:,t) = mld_tmp;
end

end