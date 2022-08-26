% Import CHL
function Chlor_a = import_CHL_OSU(lat,lon,month,ocean_mask)

load('Data/CHL.mat','CHL');

% Convert longitude
CHL.lon = convert_lon(CHL.lon);

% Reorder according to longitude
[CHL.lon,lonidx] = sort(CHL.lon);
CHL.chl = CHL.chl(:,lonidx,:);

% Format latitude and longitude
CHL.latitude = repmat(CHL.lat',1,size(CHL.chl,2),size(CHL.chl,3));
CHL.longitude = repmat(CHL.lon,size(CHL.chl,1),1,size(CHL.chl,3));

% Cut down dataset to limits of LME
lonidx = CHL.lon >= min(lon) & CHL.lon <= max(lon);
latidx = CHL.lat >= min(lat) & CHL.lat <= max(lat);
CHL.chl = double(CHL.chl(latidx,lonidx,:));
CHL.latitude = CHL.latitude(latidx,lonidx,:);
CHL.longitude = CHL.longitude(latidx,lonidx,:);

% % Create climatology
% for m = 1:12
%     chl_clim(:,:,m) = mean(CHL.chl(:,:,m:12:end),3,'omitnan');
% end

% Pre-allocate
Chlor_a = nan(length(lon),length(lat),length(month));

% Interpolate onto quarter degree grid
for t = 1:length(month)
    interp = griddedInterpolant(flipud(CHL.longitude(:,:,t))',...
        flipud(CHL.latitude(:,:,t))',flipud(CHL.chl(:,:,t))');
    lon_tmp = repmat(lon,1,length(lat));
    lat_tmp = repmat(lat',length(lon),1);
    Chlor_a(:,:,t) = interp(lon_tmp,lat_tmp);
end

% Interpolate over some gaps in CHL dataset (linear, 1-D, time), then
% remaining gaps at either end (nearest, 1-D, time)
for g = 1:length(lon)
    for h = 1:length(lat)
        if sum(~isnan(Chlor_a(g,h,:))) >= 100 % check for "too many" NaNs
            % linear interpolation
            chl_tmp = squeeze(Chlor_a(g,h,:));
            idx = ~isnan(chl_tmp);
            Chlfit = interp1(month(idx),chl_tmp(idx),month,'linear');
            % nearest neighbor interpolation
            chl_tmp = Chlfit;
            idx = ~isnan(chl_tmp);
            Chlfit = interp1(month(idx),chl_tmp(idx),month,'nearest','extrap');
            Chlor_a(g,h,:) = Chlfit;
        else
            Chlor_a(g,h,:) = NaN;
        end
    end
end

end




