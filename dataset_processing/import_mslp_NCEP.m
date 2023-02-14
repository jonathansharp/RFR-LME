% Import NCEP atmospheric pressure
function mslp = import_Pres_NCEP(lat,lon,month,ocean_mask,path)

load([path '/Data/global/NCEP.mat'],'NCEP');

% Format latitude and longitude
NCEP.latitude = repmat(NCEP.lat',size(NCEP.mslp,1),1);
NCEP.longitude = repmat(NCEP.lon,1,size(NCEP.mslp,2));

% Pre-allocate
mslp = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid to match limits of LME
for t = 1:length(month)
    interp = griddedInterpolant(fliplr(NCEP.longitude),...
        fliplr(NCEP.latitude),fliplr(NCEP.mslp(:,:,t)));
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    mslp(:,:,t) = interp(lon_tmp,lat_tmp);
end

% Convert Pascals to Atmospheres
mslp = mslp./101325;

% Remove values outside of ocean mask
for t = 1:length(month)
    mslp_tmp = mslp(:,:,t);
    mslp_tmp(~ocean_mask) = NaN;
    mslp(:,:,t) = mslp_tmp;
end

end