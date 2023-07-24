% Grid SOCAT data
% 
% This script grids observations of fCO2 and ancillary variables from five
% US large Marine Ecosystems (Alaska, California Current, Insular Pacific /
% Hawaii, Gulf of Mexico / Caribbean, and US East Coast) assembled via
% extractions from SOCATv2022 defined by latitude and longitude bounds.
% 
% Written by J.D. Sharp: 7/26/22
% Last updated by J.D. Sharp: 11/15/22
% 

%% load SOCAT structure
if ~exist('SOCAT','var')
    load('Data/socat_structure_2023','SOCAT');
end

%% display status
disp('Gridding SOCAT observations');

%% Establish latitude and longitude minima and maxima
SOCAT_grid.lim.latmin = round(min(SOCAT.latitude),0);
SOCAT_grid.lim.latmax = round(max(SOCAT.latitude),0);
SOCAT_grid.lim.lonmin = round(min(SOCAT.longitude),0);
SOCAT_grid.lim.lonmax = round(max(SOCAT.longitude),0);
SOCAT_grid.lim.monthmin = 1;
SOCAT_grid.lim.monthmax = 288;

%% remove values later than 2021
idx_21 = SOCAT.month_since_1998 > 288;
vars = fieldnames(SOCAT);
for v = 1:length(vars)
    SOCAT.(vars{v})(idx_21) = [];
end

%% Create 0.25 x 0.25 degree monthly grid
SOCAT_grid.lon = [SOCAT_grid.lim.lonmin+0.125:0.25:SOCAT_grid.lim.lonmax]';
SOCAT_grid.dim.x = length(SOCAT_grid.lon);
SOCAT_grid.lat = [SOCAT_grid.lim.latmin+0.125:0.25:SOCAT_grid.lim.latmax]';
SOCAT_grid.dim.y = length(SOCAT_grid.lat);
SOCAT_grid.month = [SOCAT_grid.lim.monthmin-0.5:1:SOCAT_grid.lim.monthmax-0.5]';
SOCAT_grid.dim.z = length(SOCAT_grid.month);

%% Add time variables
SOCAT_grid.year = repelem(1998:2022,12)';
SOCAT_grid.month_of_year = repmat(1:12,1,24)';

%% Determine bin number of each data point
[~,~,Xnum] = histcounts(SOCAT.longitude,...
    SOCAT_grid.lim.lonmin:0.25:SOCAT_grid.lim.lonmax);
[~,~,Ynum] = histcounts(SOCAT.latitude,...
    SOCAT_grid.lim.latmin:0.25:SOCAT_grid.lim.latmax);
[~,~,Znum] = histcounts(SOCAT.month_since_1998,...
    SOCAT_grid.lim.monthmin-1:SOCAT_grid.lim.monthmax);

%% Accumulate 3D grid by applying function to SOCAT values with bin numbers that match grid cells
subs = [Xnum, Ynum, Znum];
sz = [SOCAT_grid.dim.x,...
      SOCAT_grid.dim.y,...
      SOCAT_grid.dim.z];
% cruises
SOCAT_grid.count_ncruise = accumarray(subs, SOCAT.cruise, sz, @(x) numel(unique(x)), NaN);
% fCO2
SOCAT_grid.fco2_count_nobs = accumarray(subs, SOCAT.fCO2, sz, @numel, NaN);
SOCAT_grid.fco2_ave_unwtd = accumarray(subs, SOCAT.fCO2, sz, @nanmean, NaN);
SOCAT_grid.fco2_std_unwtd = accumarray(subs, SOCAT.fCO2, sz, @nanstd, NaN);
SOCAT_grid.fco2_max_unwtd = accumarray(subs, SOCAT.fCO2, sz, @max, NaN);
SOCAT_grid.fco2_min_unwtd = accumarray(subs, SOCAT.fCO2, sz, @min, NaN);
% temperature
SOCAT_grid.sst_count_nobs = accumarray(subs, SOCAT.temperature, sz, @numel, NaN);
SOCAT_grid.sst_ave_unwtd = accumarray(subs, SOCAT.temperature, sz, @nanmean, NaN);
SOCAT_grid.sst_std_unwtd = accumarray(subs, SOCAT.temperature, sz, @nanstd, NaN);
SOCAT_grid.sst_max_unwtd = accumarray(subs, SOCAT.temperature, sz, @max, NaN);
SOCAT_grid.sst_min_unwtd = accumarray(subs, SOCAT.temperature, sz, @min, NaN);
% salinity
SOCAT_grid.sss_count_nobs = accumarray(subs, SOCAT.salinity, sz, @numel, NaN);
SOCAT_grid.sss_ave_unwtd = accumarray(subs, SOCAT.salinity, sz, @nanmean, NaN);
SOCAT_grid.sss_std_unwtd = accumarray(subs, SOCAT.salinity, sz, @nanstd, NaN);
SOCAT_grid.sss_max_unwtd = accumarray(subs, SOCAT.salinity, sz, @max, NaN);
SOCAT_grid.sss_min_unwtd = accumarray(subs, SOCAT.salinity, sz, @min, NaN);

clear subs sz Xnum Ynum Znum

%% Determine cruise-weighted means and standard deviations
% If more than one cruise is represented in a given grid cell, replace
% the unweighted value with a cruise-weighted value

% Pre-allocate with unweighted values
SOCAT_grid.fco2_ave_wtd  = SOCAT_grid.fco2_ave_unwtd;
SOCAT_grid.fco2_std_wtd  = SOCAT_grid.fco2_std_unwtd;
SOCAT_grid.sst_ave_wtd  = SOCAT_grid.sst_ave_unwtd;
SOCAT_grid.sst_std_wtd  = SOCAT_grid.sst_std_unwtd;
SOCAT_grid.sss_ave_wtd  = SOCAT_grid.sss_ave_unwtd;
SOCAT_grid.sss_std_wtd  = SOCAT_grid.sss_std_unwtd;
SOCAT_grid.fco2_grid_uncert  = nan(size(SOCAT_grid.fco2_std_unwtd));
% Determine monthly means for individual years
for a = 1:SOCAT_grid.dim.x
    for b = 1:SOCAT_grid.dim.y
        for c = 1:SOCAT_grid.dim.z
            if SOCAT_grid.count_ncruise(a,b,c) > 1
                % index to specific grid cell
                idx = SOCAT.longitude >= SOCAT_grid.lon(a) - 0.125 & ...
                      SOCAT.longitude < SOCAT_grid.lon(a) + 0.125 & ...
                      SOCAT.latitude  >= SOCAT_grid.lat(b) - 0.125 & ...
                      SOCAT.latitude  < SOCAT_grid.lat(b) + 0.125 & ...
                      SOCAT.month  >= SOCAT_grid.month(c) & ...
                      SOCAT.month  < SOCAT_grid.month(c) + 1;
                cruises = unique(SOCAT.expocode(idx));
                cruiselist = SOCAT.expocode(idx);
                fco2 = SOCAT.fCO2(idx);
                temperature = SOCAT.temperature(idx);
                salinity = SOCAT.salinity(idx);
                fco22 = nan(numel(cruises),1);
                std22 = nan(numel(cruises),1);
                temperature2 = nan(numel(cruises),1);
                salinity2 = nan(numel(cruises),1);
                for k=1:numel(cruises)
                    cruiseidx = strcmp(cruiselist,cruises(k));
                    fco22(k) = mean(fco2(cruiseidx),'omitnan');
                    if length(fco2(cruiseidx)) > 1
                        std22(k) = std(fco2(cruiseidx));
                    else
                        std22(k) = NaN;
                    end
                    temperature2(k) = mean(temperature(cruiseidx),'omitnan');
                    salinity2(k) = mean(salinity(cruiseidx),'omitnan');
                end
                SOCAT_grid.fco2_ave_wtd(a,b,c)  = mean(fco22,'omitnan');
                SOCAT_grid.fco2_std_wtd(a,b,c)  = mean(std22,'omitnan');
                if any(startsWith(cruises,'3164'))
                    % checks for grid cells that include moorings; these
                    % are excluded from uncertainty estimates due to
                    % the high standard deviation among mooring observations
                    % keyboard
                else
                    SOCAT_grid.fco2_grid_uncert(a,b,c)  = ...
                        SOCAT_grid.fco2_std_unwtd(a,b,c);
                end
                SOCAT_grid.sst_ave_wtd(a,b,c)  = mean(temperature2,'omitnan');
                SOCAT_grid.sst_std_wtd(a,b,c)  = std(temperature2,[],'omitnan');
                SOCAT_grid.sss_ave_wtd(a,b,c)  = mean(salinity2,'omitnan');
                SOCAT_grid.sss_std_wtd(a,b,c)  = std(salinity2,[],'omitnan');
            end
        end
    end
end

clear a b c idx cruises cruiselist fco2 temperature salinity
clear fco22 std22 temperature2 salinity2 k cruiseidx SOCAT

%% Count number of months with at least one observation in each grid cell
% Total
SOCAT_grid.num_months = sum(~isnan(SOCAT_grid.fco2_ave_unwtd),3);
% Climatological
SOCAT_grid.fco2_ave_unwtd_clim = nan(SOCAT_grid.dim.x,SOCAT_grid.dim.y,12);
for m = 1:12
    SOCAT_grid.fco2_ave_unwtd_clim(:,:,m) = ...
        mean(SOCAT_grid.fco2_ave_unwtd(:,:,m:12:end),3,'omitnan');
end
SOCAT_grid.num_months_clim = sum(~isnan(SOCAT_grid.fco2_ave_unwtd_clim),3);

%% Determine area of each grid cell
SOCAT_grid.area_km2 = ...
(((repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1) + 0.125) - ...
    (repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1) - 0.125)) .* 110.574) .* ... % latitude distance
(((repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y) + 0.125) - ...
    (repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y) - 0.125)) .* ...
    111.320.*cosd(repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1))); % longitude distance

%% Determine sea fraction of each grid cell
load('data_to_use/ETOPO2.mat','ETOPO2');
ETOPO2.lon = convert_lon(ETOPO2.lon);
%     % limit to LME in question
%     lonidx = ETOPO2.lon >= SOCAT_grid.lim.lonmin - 360 & ...
%         ETOPO2.lon <= SOCAT_grid.lim.lonmax - 360;
%     latidx = ETOPO2.lat >= SOCAT_grid.lim.latmin & ...
%         ETOPO2.lat <= SOCAT_grid.lim.latmax;
%     ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx);
%     ETOPO2.lon = ETOPO2.lon(lonidx);
%     ETOPO2.lat = ETOPO2.lat(latidx);
% define points as land (0) or sea (1)
ETOPO2.sea = ETOPO2.bottomdepth > 0;
% determine percentage sea in each grid cell
SOCAT_grid.percent_sea = nan(size(SOCAT_grid.lat));
for a = 1:length(SOCAT_grid.lon)
    for b = 1:length(SOCAT_grid.lat)
        lonidx = find(ETOPO2.lon >= (SOCAT_grid.lon(a))-0.125 & ...
                  ETOPO2.lon < (SOCAT_grid.lon(a))+0.125);
        latidx = find(ETOPO2.lat >= SOCAT_grid.lat(b)-0.125 & ...
                  ETOPO2.lat < SOCAT_grid.lat(b)+0.125);
        SOCAT_grid.percent_sea(a,b) = ...
            sum(sum(ETOPO2.sea(lonidx,latidx)))./...
            (size(ETOPO2.sea(lonidx,latidx),1)*size(ETOPO2.sea(lonidx,latidx),2));
    end
end
clear ETOPO2 lonidx latidx a b

%% Plot the percentage of grid cells with data
figure('visible','off');
worldmap([SOCAT_grid.lim.latmin ...
    SOCAT_grid.lim.latmax],...
   [SOCAT_grid.lim.lonmin ...
    SOCAT_grid.lim.lonmax]);
setm(gca,'MapProjection','miller');
set(gca,'fontsize',16);
pcolorm(repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1),...
        repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y),...
        SOCAT_grid.num_months);
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
mycolormap = jet(16);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-1 31]);
c.TickLength = 0;
c.Label.String = 'Number of Months Represented';
cbarrow('up');
mlabel off;
plabel off;

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/obs.png');
close

% clean up
clear land c mycolormap

%% Detrend gridded pCO2 using domain mean
% Define area weights
area_weights = SOCAT_grid.area_km2;
area_weights = repmat(area_weights,1,1,size(SOCAT_grid.fco2_ave_wtd,3));
area_weights(isnan(SOCAT_grid.fco2_ave_wtd)) = NaN;
% Calculate area-weighted domain mean
SOCAT_grid.fco2_dom_mean = ...
    squeeze(sum(sum(SOCAT_grid.fco2_ave_wtd.*...
    area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
SOCAT_grid.fco2_dom_mean(SOCAT_grid.fco2_dom_mean == 0) = NaN;
% Fit trend to area weighted domain mean
[yf,yr,x] = leastsq2(SOCAT_grid.month,...
    SOCAT_grid.fco2_dom_mean,0,0,0);
% Remove difference from mean for each month
for m = 1:length(SOCAT_grid.month)
    SOCAT_grid.fco2_ave_wtd_detrend(:,:,m) = ...
        SOCAT_grid.fco2_ave_wtd(:,:,m) + ...
        (mean(yf,'omitnan') - yf(m,:));
end
% Calculate area-weighted detrended domain mean
SOCAT_grid.fco2_dom_mean_detrend = ...
    squeeze(sum(sum(SOCAT_grid.fco2_ave_wtd_detrend.*...
    area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
SOCAT_grid.fco2_dom_mean_detrend(SOCAT_grid.fco2_dom_mean_detrend == 0) = NaN;

% clean up
clear area_weights yf yr x m

%% Plot detrended gridded mean pCO2
figure('visible','off');
worldmap([SOCAT_grid.lim.latmin ...
    SOCAT_grid.lim.latmax],...
   [SOCAT_grid.lim.lonmin ...
    SOCAT_grid.lim.lonmax]);
setm(gca,'MapProjection','miller');
set(gca,'fontsize',16);
pcolorm(repmat(SOCAT_grid.lat',SOCAT_grid.dim.x,1),...
        repmat(SOCAT_grid.lon,1,SOCAT_grid.dim.y),...
        mean(SOCAT_grid.fco2_ave_wtd_detrend,3,'omitnan'));
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(parula(20));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Surface {\itf}CO_{2} (\muatm)';
cbarrow;
mlabel off;
plabel off;

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/fCO2.png');
close

% clean up
clear land c mycolormap

%% Save gridded pco2 data
if ~isfolder('Data'); mkdir('Data'); end
save('Data/socat_gridded_2023','SOCAT_grid','-v7.3');

%% clean up
clear SOCAT_grid

