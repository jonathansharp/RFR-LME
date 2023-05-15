% convert longitude from either 0:360 to -180:180 or vice-versa

function lon = convert_lon(lon)

if any(round(lon,6) > 180)
    % convert from 0:360 to -180:180
    lon(lon>180) = lon(lon>180)-360;
elseif any(round(lon,6) < 0)
    % convert from -180:180 to 0:360
    lon(lon<0) = lon(lon<0)+360;
else
    % do nothing
end