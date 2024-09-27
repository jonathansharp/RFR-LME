% Extract regional SOCAT data
% 
% This script grids observations of fCO2 and ancillary variables from five
% US large Marine Ecosystems (Alaska, California Current, Insular Pacific /
% Hawaii, Gulf of Mexico / Caribbean, and US East Coast) assembled via
% extractions from SOCATv2024 defined by latitude and longitude bounds.
% 
% Written by J.D. Sharp: 7/26/22
% Last updated by J.D. Sharp: 7/5/24
% 

%% this script defines the bounds of the eleven LMEs
define_regions_eiwg

%% extract each LME from large grid
for n = 1:length(region)

    %% load SOCAT grid
    load('Data/socat_gridded_2023','SOCAT_grid');

    %% display status
    disp(['Extracting ' region{n} ' LME from SOCAT grid']);

    %% remove observations outside general LME limits
    % determine geographic indices
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X)';
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_xmin = find(abs(SOCAT_grid.lon - min(tmp_lon))==...
        min(abs(SOCAT_grid.lon - min(tmp_lon))));
    idx_xmax = find(abs(SOCAT_grid.lon - max(tmp_lon))==...
        min(abs(SOCAT_grid.lon - max(tmp_lon))));
    idx_ymin = find(abs(SOCAT_grid.lat - min(tmp_lat))==...
        min(abs(SOCAT_grid.lat - min(tmp_lat))));
    idx_ymax = find(abs(SOCAT_grid.lat - max(tmp_lat))==...
        min(abs(SOCAT_grid.lat - max(tmp_lat))));
    % pre-allocate
    vars = fields(SOCAT_grid);
    % remove gridded observations outside region
    for v = 1:length(vars)
       if size(SOCAT_grid.(vars{v}),2) == SOCAT_grid.dim.y && ...
               ~strcmp(vars{v},'idxspc')
           SOCAT_grid.(region{n}).(vars{v}) = ...
               SOCAT_grid.(vars{v})(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
           SOCAT_grid = rmfield(SOCAT_grid,vars{v});
       end
    end
    % copy other variables
    SOCAT_grid.(region{n}).lon = SOCAT_grid.lon(idx_xmin:idx_xmax);
    SOCAT_grid.(region{n}).lat = SOCAT_grid.lat(idx_ymin:idx_ymax);
    SOCAT_grid.(region{n}).month = SOCAT_grid.month;
    SOCAT_grid.(region{n}).year = SOCAT_grid.year;
    SOCAT_grid.(region{n}).month_of_year = SOCAT_grid.month_of_year;
    SOCAT_grid.(region{n}).lim.latmin = min(SOCAT_grid.(region{n}).lat)-0.25;
    SOCAT_grid.(region{n}).lim.latmax = max(SOCAT_grid.(region{n}).lat)+0.25;
    SOCAT_grid.(region{n}).lim.lonmin = min(SOCAT_grid.(region{n}).lon)-0.25;
    SOCAT_grid.(region{n}).lim.lonmax = max(SOCAT_grid.(region{n}).lon)+0.25;
    SOCAT_grid.(region{n}).lim.monthmin = SOCAT_grid.lim.monthmin;
    SOCAT_grid.(region{n}).lim.monthmax = SOCAT_grid.lim.monthmax;
    SOCAT_grid.(region{n}).dim.x = length(SOCAT_grid.(region{n}).lon);
    SOCAT_grid.(region{n}).dim.y = length(SOCAT_grid.(region{n}).lat);
    SOCAT_grid.(region{n}).dim.z = length(SOCAT_grid.(region{n}).month);
    SOCAT_grid = rmfield(SOCAT_grid,{'lim' 'lon' 'dim' 'lat' 'month' ...
        'year' 'month_of_year' 'fco2_dom_mean' 'fco2_dom_mean_detrend'});

    % Save gridded pco2 data
    if ~isfolder(['Data/' region{n}]); mkdir(['Data/' region{n}]); end
    save(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid','-v7.3');

    % clean up
    clear tmp_lon tmp_lat idx_xmin idx_xmax idx_ymin idx_ymax vars v

    %% remove observations outside refined LME limits
    % determine index based on LME
    SOCAT_grid.(region{n}).idxspc = ...
        nan(SOCAT_grid.(region{n}).dim.x,SOCAT_grid.(region{n}).dim.y);
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    SOCAT_grid.(region{n}).idxspc = ...
        inpolygon(...
        repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y),...
        repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1),...
        tmp_lon,tmp_lat);
    SOCAT_grid.(region{n}).idxspc = ...
        repmat(SOCAT_grid.(region{n}).idxspc,1,1,SOCAT_grid.(region{n}).dim.z);
    % eliminate gridded data outside LME
    vars = fields(SOCAT_grid.(region{n}));
    for v = 1:length(vars)
       if size(SOCAT_grid.(region{n}).(vars{v}),2) == SOCAT_grid.(region{n}).dim.y && ...
               ~strcmp(vars{v},'idxspc')
           if size(SOCAT_grid.(region{n}).(vars{v}),3) == SOCAT_grid.(region{n}).dim.z
               SOCAT_grid.(region{n}).(vars{v})(~SOCAT_grid.(region{n}).idxspc) = NaN;
           else
               SOCAT_grid.(region{n}).(vars{v})(~SOCAT_grid.(region{n}).idxspc(:,:,1)) = NaN;
           end
       end
    end

    % calculate area of LME
    SOCAT_grid.(region{n}).area_tot_km2 = ...
        sum(SOCAT_grid.(region{n}).area_km2(SOCAT_grid.(region{n}).idxspc(:,:,1)).*...
        SOCAT_grid.(region{n}).percent_sea(SOCAT_grid.(region{n}).idxspc(:,:,1)),...
        'omitnan');
    SOCAT_grid.(region{n}).area_tot_km2

    % clean up
    clear tmp_lon vars v

    %% Detrend gridded pCO2 using domain mean
    % Define area weights
    area_weights = SOCAT_grid.(region{n}).area_km2;
    area_weights = repmat(area_weights,1,1,size(SOCAT_grid.(region{n}).fco2_ave_wtd,3));
    area_weights(isnan(SOCAT_grid.(region{n}).fco2_ave_wtd)) = NaN;
    % Calculate area-weighted domain mean
    SOCAT_grid.(region{n}).fco2_dom_mean = ...
        squeeze(sum(sum(SOCAT_grid.(region{n}).fco2_ave_wtd.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    SOCAT_grid.(region{n}).fco2_dom_mean(SOCAT_grid.(region{n}).fco2_dom_mean == 0) = NaN;
    % Fit trend to area weighted domain mean
    [yf,yr,x] = leastsq2(SOCAT_grid.(region{n}).month,...
        SOCAT_grid.(region{n}).fco2_dom_mean,0,0,0);
    % Remove difference from mean for each month
    for m = 1:length(SOCAT_grid.(region{n}).month)
        SOCAT_grid.(region{n}).fco2_ave_wtd_detrend(:,:,m) = ...
            SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m) + ...
            (mean(yf,'omitnan') - yf(m,:));
    end
    % Calculate area-weighted detrended domain mean
    SOCAT_grid.(region{n}).fco2_dom_mean_detrend = ...
        squeeze(sum(sum(SOCAT_grid.(region{n}).fco2_ave_wtd_detrend.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    SOCAT_grid.(region{n}).fco2_dom_mean_detrend(SOCAT_grid.(region{n}).fco2_dom_mean_detrend == 0) = NaN;

    % clean up
    clear area_weights yf yr x m

    % Save gridded pco2 data
    if ~isfolder(['Data/' region{n}]); mkdir(['Data/' region{n}]); end
    save(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid','-v7.3');

    % clean up
    clear SOCAT_grid

end