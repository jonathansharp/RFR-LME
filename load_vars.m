% Load Vars
% 
% This script grids loads satellite variables for five US large Marine
% Ecosystems (Alaska, California Current, Insular Pacific / Hawaii, Gulf of
% Mexico / Caribbean, and US East Coast) and regrids them in preparation
% for clustering and machine learning algorithm fits.
% 
% Written by J.D. Sharp: 8/19/22
% Last updated by J.D. Sharp: 8/23/22
% 

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 4%1:length(region)

    %% display status
    disp(['Loading predictor variables (' region{n} ')']);

    %% landmask
    lon_tmp = repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y);
    lon_tmp(lon_tmp > 180) = lon_tmp(lon_tmp > 180) - 360;
    lat_tmp = repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1);
    ocean_mask = ~landmask(lat_tmp,lon_tmp);

    %% Obtain sea surface salinity from ECCO reanalysis
    SOCAT.([region{n} '_grid']).SSS = ...
        import_SSS_ECCO(SOCAT.([region{n} '_grid']).lat,...
                        SOCAT.([region{n} '_grid']).lon,ocean_mask);
    
%     % test plot
%     figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
%     pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
%         repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
%         SOCAT.([region{n} '_grid']).SSS(:,:,1));
%     colorbar;

    %% Obtain sea surface height from ECCO reanalysis
    SOCAT.([region{n} '_grid']).SSH = ...
        import_SSH_ECCO(SOCAT.([region{n} '_grid']).lat,...
                        SOCAT.([region{n} '_grid']).lon,ocean_mask);
    
    % test plot
    figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
    pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
        repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
        SOCAT.([region{n} '_grid']).SSH(:,:,171));
    colorbar;

%% Obtain sea surface temperature from OISSTv2
    SOCAT.([region{n} '_grid']).SST = ...
        import_SST_OISST(SOCAT.([region{n} '_grid']).lat,...
                        SOCAT.([region{n} '_grid']).lon,ocean_mask);
    
%     % test plot
%     figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
%     pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
%         repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
%         SOCAT.([region{n} '_grid']).SST(:,:,19));
%     colorbar;

%% Obtain sea surface chlorophyll from satellite measurements
    SOCAT.([region{n} '_grid']).CHL = ...
        import_CHL_OSU(SOCAT.([region{n} '_grid']).lat,...
                        SOCAT.([region{n} '_grid']).lon,...
                        SOCAT.([region{n} '_grid']).month,ocean_mask);
    
%     % test plot
%     figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
%     pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
%         repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
%         log10(SOCAT.([region{n} '_grid']).CHL(:,:,223)));
%     colorbar;

end

%% Obtain wind speed from ERA5 re-analysis


%% Obtain wind speed from NCEP re-analysis
disp('Obtaining NCEP re-analysis winds');
load('Data/NCEPw.mat');
% Cut down dataset to CCS limits
lonidx = NCEPw.lon(:,1,1) >= lonmin & NCEPw.lon(:,1,1) <= lonmax;
latidx = NCEPw.lat(1,:,1)' >= latmin & NCEPw.lat(1,:,1)' <= latmax;
NCEPw.speed = NCEPw.speed(lonidx,latidx,:);
NCEPw.u = NCEPw.u(lonidx,latidx,:);
NCEPw.v = NCEPw.v(lonidx,latidx,:);
NCEPw.lat = NCEPw.lat(lonidx,latidx,:);
NCEPw.lon = NCEPw.lon(lonidx,latidx,:);
% Match time frame of SOCAT data
NCEPw.date = datevec(NCEPw.date);
NCEPw.month_since_1998 = (NCEPw.date(:,1)-1998).*12 + NCEPw.date(:,2);
idx = NCEPw.month_since_1998 >= min(SOCATv2021_grid.month_since_1998) & ...
      NCEPw.month_since_1998 <= max(SOCATv2021_grid.month_since_1998);
NCEPw.speed = NCEPw.speed(:,:,idx);
NCEPw.u = NCEPw.u(:,:,idx);
NCEPw.v = NCEPw.v(:,:,idx);
NCEPw.lat = NCEPw.lat(:,:,idx);
NCEPw.lon = NCEPw.lon(:,:,idx);
disp('Interpolating wind speed to SOCAT grid');
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(NCEPw.lon(:,:,t),NCEPw.lat(:,:,t),NCEPw.speed(:,:,t));
    SOCATv2021_grid.wind_speed_NCEP(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
    interp = griddedInterpolant(NCEPw.lon(:,:,t),NCEPw.lat(:,:,t),NCEPw.u(:,:,t));
    SOCATv2021_grid.u10_NCEP(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
    interp = griddedInterpolant(NCEPw.lon(:,:,t),NCEPw.lat(:,:,t),NCEPw.v(:,:,t));
    SOCATv2021_grid.v10_NCEP(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
end

% figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.wind_speed(:,:,1)); colorbar; caxis([0 12]);

clear t idx latidx lonidx interp NCEPw

%% Obtain wind speed from CCMP re-analysis
disp('Obtaining CCMP re-analysis winds');
load('Data/CCMP.mat');
% Cut down dataset to CCS limits
lonidx = CCMP.lon(:,1,1) >= lonmin & CCMP.lon(:,1,1) <= lonmax;
latidx = CCMP.lat(1,:,1) >= latmin & CCMP.lat(1,:,1) <= latmax;
CCMP.speed = CCMP.speed(lonidx,latidx,:);
CCMP.u10 = CCMP.u10(lonidx,latidx,:);
CCMP.v10 = CCMP.v10(lonidx,latidx,:);
CCMP.lat = CCMP.lat(lonidx,latidx,:);
CCMP.lon = CCMP.lon(lonidx,latidx,:);
% Match time frame of SOCAT data
CCMP.speed = repmat(CCMP.speed,1,1,23);
CCMP.u10 = repmat(CCMP.u10,1,1,23);
CCMP.v10 = repmat(CCMP.v10,1,1,23);
disp('Interpolating wind speed to SOCAT grid');
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(CCMP.lon,CCMP.lat,CCMP.speed(:,:,t));
    SOCATv2021_grid.wind_speed_CCMP(:,:,t) = interp(SOCATv2021_grid.lon,SOCATv2021_grid.lat);
    interp = griddedInterpolant(CCMP.lon,CCMP.lat,CCMP.u10(:,:,t));
    SOCATv2021_grid.u10_CCMP(:,:,t) = interp(SOCATv2021_grid.lon,SOCATv2021_grid.lat);
    interp = griddedInterpolant(CCMP.lon,CCMP.lat,CCMP.v10(:,:,t));
    SOCATv2021_grid.v10_CCMP(:,:,t) = interp(SOCATv2021_grid.lon,SOCATv2021_grid.lat);
end

% figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.wind_speed(:,:,1)); colorbar; caxis([0 12]);

clear t idx latidx lonidx interp CCMP

%% Obtain bathymetry from ETOPO2
disp('Obtaining ETOPO2 bathymetry');
load('Data/ETOPO2.mat');
ETOPO2.latitude = repmat(ETOPO2.lat',size(ETOPO2.bottomdepth,1),1);
ETOPO2.longitude = repmat(ETOPO2.lon,1,size(ETOPO2.bottomdepth,2));
% Cut down dataset to CCS limits
lonidx = ETOPO2.lon >= lonmin-360 & ETOPO2.lon <= lonmax-360;
latidx = ETOPO2.lat >= latmin & ETOPO2.lat <= latmax;
ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx,:);
ETOPO2.latitude = ETOPO2.latitude(lonidx,latidx,:);
ETOPO2.longitude = ETOPO2.longitude(lonidx,latidx,:);
ETOPO2.longitude = ETOPO2.longitude + 360; % Convert to 360 degrees longitude
% Interpolate onto SOCAT grid
disp('Interpolating bathymetry to SOCAT grid');
interp = griddedInterpolant(ETOPO2.longitude,ETOPO2.latitude,ETOPO2.bottomdepth);
SOCATv2021_grid.bottomdepth = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
% Extend to time frame of SOCAT data
SOCATv2021_grid.bottomdepth = repmat(SOCATv2021_grid.bottomdepth,1,1,max(size(SOCATv2021_grid.month_since_1998)));

% figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.bottomdepth(:,:,1)); colorbar;

clear path latidx lonidx interp ETOPO2

%% Obtain mixed layer depth from HYCOM model
disp('Obtaining HYCOM mixed layer depth');
load('Data/MLD.mat');
MLD.latitude = repmat(MLD.lat',1,size(MLD.mld,2),size(MLD.mld,3));
MLD.longitude = repmat(MLD.lon,size(MLD.mld,1),1,size(MLD.mld,3));
% Cut down dataset to CCS limits
lonidx = MLD.lon >= lonmin-360 & MLD.lon <= lonmax-360;
latidx = MLD.lat >= latmin & MLD.lat <= latmax;
MLD.mld = MLD.mld(latidx,lonidx,:);
MLD.latitude = MLD.latitude(latidx,lonidx,:);
MLD.longitude = MLD.longitude(latidx,lonidx,:);
MLD.longitude = MLD.longitude + 360; % Convert to 360 degrees longitude
% Match time frame of SOCAT data
MLD.date = datevec(MLD.time);
MLD.month_since_1998 = (MLD.date(:,1)-1998).*12 + MLD.date(:,2);
idx = MLD.month_since_1998 >= min(SOCATv2021_grid.month_since_1998) & ...
      MLD.month_since_1998 <= max(SOCATv2021_grid.month_since_1998);
MLD.mld = MLD.mld(:,:,idx);
MLD.latitude = MLD.latitude(:,:,idx);
MLD.longitude = MLD.longitude(:,:,idx);
% Interpolate onto SOCAT grid
disp('Interpolating mixed layer depth to SOCAT grid');
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(flipud(MLD.longitude(:,:,t))',flipud(MLD.latitude(:,:,t))',flipud(MLD.mld(:,:,t))');
    SOCATv2021_grid.MLD(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
end
% Interpolate over some gaps in MLD dataset (2-D, lat and lon)
for t = 1:max(SOCATv2021_grid.month_since_1998)
    idx = ~isnan(SOCATv2021_grid.MLD(:,:,t)) & ~SOCATv2021_grid.mask_land;
    lon_tmp = SOCATv2021_grid.longitude(:,:,t);
    lat_tmp = SOCATv2021_grid.latitude(:,:,t);
    mld_tmp = SOCATv2021_grid.MLD(:,:,t);
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),mld_tmp(idx));
    idx = isnan(SOCATv2021_grid.MLD(:,:,t)) & ~SOCATv2021_grid.mask_land;
    mld_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    SOCATv2021_grid.MLD(:,:,t) = mld_tmp;
end
% Eliminate lake values
SOCATv2021_grid.MLD(SOCATv2021_grid.latitude > 40 & ...
    SOCATv2021_grid.longitude(:,:,1) > 240) = NaN;
% % Replace values less than zero with value close to zero
% SOCATv2021_grid.MLD(SOCATv2021_grid.MLD<0) = 0.001;

% figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.MLD(:,:,1)); colorbar; caxis([0 150]);

clear latidx lonidx t idx lat_tmp lon_tmp mld_tmp interp MLD

% %% Obtain atmospheric pCO2 from Jena CarboScope
% disp('Obtaining Jena CarboScope atmospheric pCO2');
% importCarboScope
% CarboScope.latitude = repmat(CarboScope.lat',size(CarboScope.paCO2,1),1,size(CarboScope.paCO2,3));
% CarboScope.longitude = repmat(CarboScope.lon,1,size(CarboScope.paCO2,2),size(CarboScope.paCO2,3));
% % Cut down dataset to CCS limits
% lonidx = CarboScope.lon >= lonmin-360 & CarboScope.lon <= lonmax-360;
% latidx = CarboScope.lat >= latmin & CarboScope.lat <= latmax;
% CarboScope.paCO2 = CarboScope.paCO2(lonidx,latidx,:);
% CarboScope.latitude = CarboScope.latitude(lonidx,latidx,:);
% CarboScope.longitude = CarboScope.longitude(lonidx,latidx,:);
% CarboScope.longitude = CarboScope.longitude + 360; % Convert to 360 degrees longitude
% CarboScope.latitude = CarboScope.latitude(:,:,1);
% CarboScope.longitude = CarboScope.longitude(:,:,1);
% % Match time frame of SOCAT data
% CarboScope.date = datevec(CarboScope.date);
% CarboScope.month_since_1998 = (CarboScope.date(:,1)-1998).*12 + CarboScope.date(:,2);
% idx = CarboScope.month_since_1998 >= min(SOCATv2021_grid.month_since_1998) & ...
%       CarboScope.month_since_1998 <= max(SOCATv2021_grid.month_since_1998);
% CarboScope.date = CarboScope.date(idx,1:3);
% CarboScope.month_since_1998 = CarboScope.month_since_1998(idx);
% CarboScope.paCO2 = CarboScope.paCO2(:,:,idx);
% % Determine monthly averages
% CarboScope.paCO2_mon = ...
%     nan(size(CarboScope.longitude,1),size(CarboScope.latitude,2),max(SOCATv2021_grid.month_since_1998));
% for t = 1:max(SOCATv2021_grid.month_since_1998)
%     idx = CarboScope.month_since_1998 == t;
%     CarboScope.paCO2_mon(:,:,t) = mean(CarboScope.paCO2(:,:,idx),3,'omitnan');
% end
% % Interpolate onto SOCAT grid
% disp('Interpolating atmospheric pCO2 to SOCAT grid');
% for t = 1:max(SOCATv2021_grid.month_since_1998)
%     interp = griddedInterpolant(CarboScope.longitude(:,:),CarboScope.latitude(:,:),CarboScope.paCO2_mon(:,:,t));
%     SOCATv2021_grid.paCO2(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
% end
% SOCATv2021_grid.paCO2 = double(SOCATv2021_grid.paCO2);
% 
% % figure; pcolor(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.paCO2(:,:,1)); colorbar; %caxis([0 150]);
% 
% clear latidx lonidx t idx interp

%% Obtain atmospheric pressure from NCEP
disp('Obtaining NCEP atmospheric pressure');
% Read netcdf file
NCEP = netcdfreader('Data/mslp.mon.mean.nc');
NCEP.latitude = repmat(NCEP.lat',size(NCEP.mslp,1),1);
NCEP.longitude = repmat(NCEP.lon,1,size(NCEP.mslp,2));
% Calculate date
NCEP.date = datevec(datenum([repmat(1800,size(NCEP.time,1),1) ones(size(NCEP.time,1),1) ones(size(NCEP.time,1),1) ...
    NCEP.time zeros(size(NCEP.time,1),1) zeros(size(NCEP.time,1),1)]));
NCEP.months_since_1998 = (NCEP.date(:,1) - 1998).*12 + NCEP.date(:,2);
% Cut down dataset to CCS limits
lonidx = NCEP.lon >= lonmin & NCEP.lon <= lonmax;
latidx = NCEP.lat >= latmin & NCEP.lat <= latmax;
NCEP.mslp = NCEP.mslp(lonidx,latidx,:);
NCEP.latitude = NCEP.latitude(lonidx,latidx);
NCEP.longitude = NCEP.longitude(lonidx,latidx);
% Index to 1998-2019
index = NCEP.months_since_1998 > 0 & NCEP.months_since_1998 <= max(SOCATv2021_grid.month_since_1998);
NCEP.months_since_1998 = NCEP.months_since_1998(index);
NCEP.mslp = NCEP.mslp(:,:,index);
% Interpolate onto SOCAT grid
disp('Interpolating atmospheric to SOCAT grid');
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(fliplr(NCEP.longitude(:,:)),fliplr(NCEP.latitude(:,:)),fliplr(NCEP.mslp(:,:,t)));
    SOCATv2021_grid.mslp(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
end
% Convert pascals to atmospheres
SOCATv2021_grid.mslp = double(SOCATv2021_grid.mslp)./101325;

clear t interp lonidx latidx index NCEP

%% Obtain atmospheric pCO2 from NOAA MBL product
disp('Obtaining NOAA MBL atmospheric pCO2');
% Open and scan file
file = fopen('Data/MBL_1998_2020.txt');
NOAA_MBL = textscan(file,'%f','Delimiter',',','CommentStyle','#');
fclose(file);
NOAA_MBL = reshape(cell2mat(NOAA_MBL),83,[]);
% Process data within file
MBL.year = repmat(NOAA_MBL(1,:),(size(NOAA_MBL,1)-1)./2,1);
MBL.CO2 = NOAA_MBL(2:2:end,:);
MBL.err = NOAA_MBL(3:2:end,:);
% Define latitudes
latsin = [-1.00  -0.95  -0.90  -0.85  -0.80  -0.75  -0.70  -0.65  -0.60 ...
          -0.55  -0.50  -0.45  -0.40  -0.35  -0.30  -0.25  -0.20  -0.15 ...
          -0.10  -0.05   0.00   0.05   0.10   0.15   0.20   0.25   0.30 ...
           0.35   0.40   0.45   0.50   0.55   0.60   0.65   0.70   0.75 ...
           0.80   0.85   0.90   0.95   1.00];
MBL.lat = asind(latsin)';
% Define year fraction interval based on monthly slices
MBL.year_mon = repmat(min(min(MBL.year)):1/12:max(max(MBL.year)),size(MBL.CO2,1),1);
MBL.year_mon = MBL.year_mon(:,1:end-1);
MBL.year_mon = MBL.year_mon + 1/24;
% Interpolate to monthly time slices 
MBL.CO2_mon = nan(size(MBL.year_mon));
MBL.err_mon = nan(size(MBL.year_mon));
for l = 1:size(MBL.CO2,1)
    MBL.CO2_mon(l,:) = interp1(MBL.year(l,:),MBL.CO2(l,:),MBL.year_mon(l,:));
    MBL.err_mon(l,:) = interp1(MBL.year(l,:),MBL.err(l,:),MBL.year_mon(l,:));
end
% Interpolate values to monthly interval
SOCATv2021_grid.xCO2_atm = nan(size(SOCATv2021_grid.lat,2),size(MBL.CO2_mon,2));
SOCATv2021_grid.xCO2_err_atm = nan(size(SOCATv2021_grid.lat,2),size(MBL.CO2_mon,2));
for m = 1:size(MBL.CO2_mon,2)
    SOCATv2021_grid.xCO2_atm(:,m) = interp1(MBL.lat,MBL.CO2_mon(:,m),SOCATv2021_grid.lat(1,:)');
    SOCATv2021_grid.xCO2_err_atm(:,m) = interp1(MBL.lat,MBL.err_mon(:,m),SOCATv2021_grid.lat(1,:)');
end
% % Extend to 2020 using linear fit for each month
% SOCATv2021_grid.xCO2_atm = [SOCATv2021_grid.xCO2_atm nan(180,12)];
% SOCATv2021_grid.xCO2_err_atm = [SOCATv2021_grid.xCO2_err_atm nan(180,12)];
% for m = 1:12 % for each month
%     for l = 1:size(SOCATv2021_grid.xCO2_atm,1) % for each latitude
%         fit = polyfit(SOCATv2021_grid.month_since_1998(m:12:264),...
%             squeeze(SOCATv2021_grid.xCO2_atm(l,m:12:264)),1);
%         SOCATv2021_grid.xCO2_atm(l,m+264) = fit(1)*(m+264) + fit(2);
%         fit = polyfit(SOCATv2021_grid.month_since_1998(m:12:264),...
%             squeeze(SOCATv2021_grid.xCO2_err_atm(l,m:12:264)),1);
%         SOCATv2021_grid.xCO2_err_atm(l,m+264) = fit(1)*(m+264) + fit(2);
%     end
% end
% Replicate across longitudes
SOCATv2021_grid.xCO2_atm = repmat(permute(SOCATv2021_grid.xCO2_atm,[3 1 2]),size(SOCATv2021_grid.lon,1),1,1);
SOCATv2021_grid.xCO2_err_atm = repmat(permute(SOCATv2021_grid.xCO2_err_atm,[3 1 2]),size(SOCATv2021_grid.lon,1),1,1);
% Calculate pCO2 from xCO2 (with vapor pressure correction)
SOCATv2021_grid.vapor_pressure = vpress(SOCATv2021_grid.SSS,SOCATv2021_grid.SST);
SOCATv2021_grid.pCO2_atm = SOCATv2021_grid.xCO2_atm.*...
    (SOCATv2021_grid.mslp-SOCATv2021_grid.vapor_pressure);

clear ans file l t z m MBL* latsin NOAA_MBL
