% Import ERA5 Winds

load('Data/ERA5.mat');

ERA5.latitude = repmat(ERA5.lat',size(ERA5.speed,1),1,size(ERA5.speed,3));
ERA5.longitude = repmat(ERA5.lon+360,1,size(ERA5.speed,2),size(ERA5.speed,3));
% Match time frame of SOCAT data
ERA5.date = datevec(ERA5.date);
ERA5.month_since_1998 = (ERA5.date(:,1)-1998).*12 + ERA5.date(:,2);
idx = ERA5.month_since_1998 >= min(SOCATv2021_grid.month_since_1998) & ...
      ERA5.month_since_1998 <= max(SOCATv2021_grid.month_since_1998);
ERA5.speed = ERA5.speed(:,:,idx);
ERA5.u10 = ERA5.u10(:,:,idx);
ERA5.v10 = ERA5.v10(:,:,idx);
ERA5.latitude = ERA5.latitude(:,:,idx);
ERA5.longitude = ERA5.longitude(:,:,idx);
disp('Interpolating wind speed to SOCAT grid');
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.speed(:,:,t)));
    SOCATv2021_grid.wind_speed(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.u10(:,:,t)));
    SOCATv2021_grid.u10(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.v10(:,:,t)));
    SOCATv2021_grid.v10(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
end

% figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.wind_speed(:,:,1)); colorbar; caxis([0 12]);

clear t idx interp ERA5