% Grid SOCAT data but exclude mooring observations
% 
% This script grids observations of fCO2 and ancillary variables into
% eleven US large Marine Ecosystems assembled via extractions from
% SOCATv2024 defined by latitude and longitude bounds.
% 
% Written by J.D. Sharp: 7/26/23
% Last updated by J.D. Sharp: 10/17/24
% 

%% load SOCAT structure
if ~exist('SOCAT','var')
    load(['Data/socat_structure_' num2str(rfr_lme_year)],'SOCAT');
end

%% display status
disp('Gridding SOCAT observations (excluding moorings)');

%% Establish latitude and longitude minima and maxima
SOCAT_grid.lim.latmin = round(min(SOCAT.latitude),0);
SOCAT_grid.lim.latmax = round(max(SOCAT.latitude),0);
SOCAT_grid.lim.lonmin = round(min(SOCAT.longitude),0);
SOCAT_grid.lim.lonmax = round(max(SOCAT.longitude),0);
SOCAT_grid.lim.monthmin = 1;
SOCAT_grid.lim.monthmax = 312;

%% remove values later than 2023
idx_23 = SOCAT.month_since_1998 > (rfr_lme_year-1998)*12;
vars = fieldnames(SOCAT);
for v = 1:length(vars)
    SOCAT.(vars{v})(idx_23) = [];
end
clear idx_23

%% Create 0.25 x 0.25 degree monthly grid
SOCAT_grid.lon = [SOCAT_grid.lim.lonmin+0.125:0.25:SOCAT_grid.lim.lonmax]';
SOCAT_grid.dim.x = length(SOCAT_grid.lon);
SOCAT_grid.lat = [SOCAT_grid.lim.latmin+0.125:0.25:SOCAT_grid.lim.latmax]';
SOCAT_grid.dim.y = length(SOCAT_grid.lat);
SOCAT_grid.month = [SOCAT_grid.lim.monthmin-0.5:1:SOCAT_grid.lim.monthmax-0.5]';
SOCAT_grid.dim.z = length(SOCAT_grid.month);

%% Add time variables
SOCAT_grid.year = repelem(1998:rfr_lme_year-1,12)';
SOCAT_grid.month_of_year = repmat(1:12,1,rfr_lme_year-1998)';

%% exclude mooring observations
plat = unique(SOCAT.expocode);
moor_idx = ~cellfun(@isempty,regexp(plat,'3164')); % this determines which expocodes are moorings
moor_idx = moor_idx | ~cellfun(@isempty,regexp(plat,'187F')); % adds something else...
test_plat = plat(moor_idx);
idx = ~ismember(SOCAT.expocode,test_plat);
vars = fieldnames(SOCAT);
for v = 1:length(vars)
    SOCAT.(vars{v}) = SOCAT.(vars{v})(idx);
end
clear plat moor_idx test_plat idx vars v

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
% Obtain bathymetry from ETOPOv2022
ETOPO.lon = ncread('data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc','lon');
ETOPO.lat = ncread('data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc','lat');
ETOPO.bottomdepth = ncread('data_to_use/ETOPO_2022_v1_60s_N90W180_bed.nc','z');
% Convert longitude
ETOPO.lon = convert_lon(ETOPO.lon);
% define points as land (0) or sea (1)
ETOPO.sea = ETOPO.bottomdepth < 0;
% determine percentage sea in each grid cell
SOCAT_grid.percent_sea = nan(size(SOCAT_grid.lat));
for a = 1:length(SOCAT_grid.lon)
    for b = 1:length(SOCAT_grid.lat)
        lonidx = find(ETOPO.lon >= (SOCAT_grid.lon(a))-0.125 & ...
                  ETOPO.lon < (SOCAT_grid.lon(a))+0.125);
        latidx = find(ETOPO.lat >= SOCAT_grid.lat(b)-0.125 & ...
                  ETOPO.lat < SOCAT_grid.lat(b)+0.125);
        SOCAT_grid.percent_sea(a,b) = ...
            sum(sum(ETOPO.sea(lonidx,latidx)))./...
            (size(ETOPO.sea(lonidx,latidx),1)*size(ETOPO.sea(lonidx,latidx),2));
    end
end
clear ETOPO lonidx latidx a b

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

%% Save gridded pco2 data
if ~isfolder('Data'); mkdir('Data'); end
save(['Data/socat_gridded_' num2str(rfr_lme_year-1) '_no_moorings'],'SOCAT_grid','-v7.3');

%% clean up
clear SOCAT_grid
