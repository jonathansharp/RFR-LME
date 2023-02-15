% Import CHL (NASA)
function Chlor_a = import_CHL_NASA(lat,lon,month,ocean_mask,path)

load([path '/data_to_use/CHL_NASA.mat'],'CHL');

% Convert longitude
CHL.lon = convert_lon(CHL.lon);

% Reorder according to longitude
[CHL.lon,lonidx] = sort(CHL.lon);
CHL.chl = CHL.chl(lonidx,:,:);

% Format latitude and longitude
CHL.latitude = repmat(CHL.lat',size(CHL.chl,1),1);
CHL.longitude = repmat(CHL.lon,1,size(CHL.chl,2));

% Cut down dataset to limits of LME
lonidx = CHL.lon >= min(lon) & CHL.lon <= max(lon);
latidx = CHL.lat >= min(lat) & CHL.lat <= max(lat);
CHL.chl = double(CHL.chl(lonidx,latidx,:));
CHL.latitude = repmat(CHL.latitude(lonidx,latidx),1,1,length(CHL.time));
CHL.longitude = repmat(CHL.longitude(lonidx,latidx),1,1,length(CHL.time));

% Pre-allocate
Chlor_a = nan(length(lon),length(lat),length(month));

% Interpolate over some gaps in CHL dataset (linear, 1-D, time), then
% remaining gaps at either end (nearest, 1-D, time)
for g = 1:size(CHL.longitude,1)
    for h = 1:size(CHL.latitude,2)
        if sum(~isnan(CHL.chl(g,h,:))) >= 100 % check for "too many" NaNs
            % linear interpolation
            chl_tmp = squeeze(CHL.chl(g,h,:));
            idx = ~isnan(chl_tmp);
            Chlfit = interp1(month(idx),chl_tmp(idx),month,'linear');
            % nearest neighbor interpolation
            chl_tmp = Chlfit;
            idx = ~isnan(chl_tmp);
            Chlfit = interp1(month(idx),chl_tmp(idx),month,'nearest','extrap');
            CHL.chl(g,h,:) = Chlfit;
        else
            CHL.chl(g,h,:) = NaN;
        end
    end
end

% Interpolate onto quarter degree grid
for t = 1:length(month)
    % Index where CHL is both true
    idx = ~isnan(CHL.chl(:,:,t));
    % Get teporary lat, lon, CHL
    lon_tmp = CHL.longitude(:,:,t);
    lat_tmp = CHL.latitude(:,:,t);
    chl_tmp = CHL.chl(:,:,t);
    % Create interpolant over than range
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),chl_tmp(idx));
    % Get teporary lat and lon
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    % interpolate to grid
    Chlor_a(:,:,t) = interp(lon_tmp,lat_tmp);
end

% Remove values outside of ocean mask
for t = 1:length(month)
    chl_tmp = Chlor_a(:,:,t);
    chl_tmp(~ocean_mask) = NaN;
    Chlor_a(:,:,t) = chl_tmp;
end

end
