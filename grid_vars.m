% Grid LMEs
% 
% This script grids observations of fCO2 and ancillary variables from five
% US large Marine Ecosystems (Alaska, California Current, Insular Pacific /
% Hawaii, Gulf of Mexico / Caribbean, and US East Coast) assembled via
% extractions from SOCATv2022 defined by latitude and longitude bounds.
% 
% Written by J.D. Sharp: 7/26/22
% Last updated by J.D. Sharp: 8/23/22

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 1:length(region)

    %% display status
    disp(['Downloading SOCAT Data (' region{n} ')']);
    
    %% load regional SOCAT file
    load(['SOCATv2022/SOCATv2022_' region{n} '.mat']);
    
    %% assemble regional variables into a structure
    SOCAT_struct;
    
    %% remove observations before 1998
    idxyr = SOCAT.(region{n}).year >= 1998;
    vars = fieldnames(SOCAT.(region{n}));
    for v = 1:numel(vars)
            tempvar = SOCAT.(region{n}).(string(vars(v)));
            SOCAT.(region{n}).(string(vars(v))) = tempvar(idxyr);
    end
    clear v tempvar idxyr vars
    
    %% determine months since Jan 1 1998
    SOCAT.(region{n}).month_since_1998 = ...
        (SOCAT.(region{n}).year-1998).*12 + SOCAT.(region{n}).month;
    
    %% remove observations with flags other than 2 and A/B/C/D
    idxflag = SOCAT.(region{n}).fCO2_flag == 2 & ...
        (strcmp(SOCAT.(region{n}).flag,'A') | ...
         strcmp(SOCAT.(region{n}).flag,'B') | ...
         strcmp(SOCAT.(region{n}).flag,'C') | ...
         strcmp(SOCAT.(region{n}).flag,'D'));
    vars = fieldnames(SOCAT.(region{n}));
    for v = 1:numel(vars)
            tempvar = SOCAT.(region{n}).(char(vars(n)));
            SOCAT.(region{n}).(char(vars(n))) = tempvar(idxflag);
    end
    clear v idxflag vars tempvar
    
    %% obtain unique integers for each expocode
    SOCAT.(region{n}).cruise = ...
        nan(size(SOCAT.(region{n}).expocode));
    cruiselist = unique(SOCAT.(region{n}).expocode);
    for c = 1:numel(cruiselist)
        idx = strcmp(cruiselist(c),SOCAT.(region{n}).expocode);
        SOCAT.(region{n}).cruise(idx) = c;
    end
    clear c cruiselist idx
    
    %% visualize number of observations
    time = datetime(SOCAT.(region{n}).year,...
        SOCAT.(region{n}).month,SOCAT.(region{n}).day);
    figure; histogram(time);
    ylabel('Number of observations');
    xlabel('Year');
    clear time

    % save figure
    exportgraphics(gcf,['Figures/' region{n} '_hist.png']);
    close
    
    %% display status
    disp(['Gridding SOCAT observations (' region{n} ')']);
    
    %% Establish latitude and longitude minima and maxima
    SOCAT_grid.(region{n}).lim.latmin = round(min(SOCAT.(region{n}).latitude),0);
    SOCAT_grid.(region{n}).lim.latmax = round(max(SOCAT.(region{n}).latitude),0);
    SOCAT_grid.(region{n}).lim.lonmin = round(min(SOCAT.(region{n}).longitude),0);
    SOCAT_grid.(region{n}).lim.lonmax = round(max(SOCAT.(region{n}).longitude),0);
    SOCAT_grid.(region{n}).lim.monthmin = 1;
    SOCAT_grid.(region{n}).lim.monthmax = 288;
    
    %% Create 0.25 x 0.25 degree grid
%     [SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,SOCAT_grid.(region{n}).month] = ...
%         meshgrid(lim.latmin+0.125:0.25:lim.latmax,lim.lonmin+0.125:0.25:lim.lonmax,lim.monthmin-0.5:1:lim.monthmax-0.5);
%     dim.x = size(SOCAT_grid.(region{n}).lon,1);
%     dim.y = size(SOCAT_grid.(region{n}).lat,2);
%     dim.z = size(SOCAT_grid.(region{n}).month,3);

    SOCAT_grid.(region{n}).lon = [SOCAT_grid.(region{n}).lim.lonmin+0.125:0.25:SOCAT_grid.(region{n}).lim.lonmax]';
    SOCAT_grid.(region{n}).dim.x = length(SOCAT_grid.(region{n}).lon);
    SOCAT_grid.(region{n}).lat = [SOCAT_grid.(region{n}).lim.latmin+0.125:0.25:SOCAT_grid.(region{n}).lim.latmax]';
    SOCAT_grid.(region{n}).dim.y = length(SOCAT_grid.(region{n}).lat);
    SOCAT_grid.(region{n}).month = [SOCAT_grid.(region{n}).lim.monthmin-0.5:1:SOCAT_grid.(region{n}).lim.monthmax-0.5]';
    SOCAT_grid.(region{n}).dim.z = length(SOCAT_grid.(region{n}).month);

    %% Add time variables
    % year
    SOCAT_grid.(region{n}).year = repelem(1998:2021,12)';
    % month of year
    SOCAT_grid.(region{n}).month_of_year = repmat(1:12,1,24)';
    % date
%     SOCAT_grid.(region{n}).date = ...
%         repmat(permute(...
%         datenum([squeeze(SOCAT_grid.(region{n}).year(1,1,:)) ...
%         squeeze(SOCAT_grid.(region{n}).month(1,1,:)+0.5) ...
%         repmat(15,SOCAT_grid.(region{n}).dim.z,1)]),...
%         [3 2 1]),SOCAT_grid.(region{n}).dim.x,...
%         SOCAT_grid.(region{n}).dim.y,1);

    %% Determine bin number of each data point
    [~,~,Xnum] = histcounts(SOCAT.(region{n}).longitude,...
        SOCAT_grid.(region{n}).lim.lonmin:0.25:SOCAT_grid.(region{n}).lim.lonmax);
    [~,~,Ynum] = histcounts(SOCAT.(region{n}).latitude,...
        SOCAT_grid.(region{n}).lim.latmin:0.25:SOCAT_grid.(region{n}).lim.latmax);
    [~,~,Znum] = histcounts(SOCAT.(region{n}).month_since_1998,...
        SOCAT_grid.(region{n}).lim.monthmin-1:SOCAT_grid.(region{n}).lim.monthmax);
    
    %% Accumulate 3D grid by applying function to SOCAT values with bin numbers that match grid cells
    subs = [Xnum, Ynum, Znum];
    sz = [SOCAT_grid.(region{n}).dim.x,...
          SOCAT_grid.(region{n}).dim.y,...
          SOCAT_grid.(region{n}).dim.z];
    % cruises
    SOCAT_grid.(region{n}).count_ncruise = accumarray(subs, SOCAT.(region{n}).cruise, sz, @(x) numel(unique(x)), NaN);
    % fCO2
    SOCAT_grid.(region{n}).fco2_count_nobs = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @numel, NaN);
    SOCAT_grid.(region{n}).fco2_ave_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @nanmean, NaN);
    SOCAT_grid.(region{n}).fco2_std_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @nanstd, NaN);
    SOCAT_grid.(region{n}).fco2_max_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @max, NaN);
    SOCAT_grid.(region{n}).fco2_min_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @min, NaN);
    % temperature
    SOCAT_grid.(region{n}).sst_count_nobs = accumarray(subs, SOCAT.(region{n}).temperature, sz, @numel, NaN);
    SOCAT_grid.(region{n}).sst_ave_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @nanmean, NaN);
    SOCAT_grid.(region{n}).sst_std_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @nanstd, NaN);
    SOCAT_grid.(region{n}).sst_max_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @max, NaN);
    SOCAT_grid.(region{n}).sst_min_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @min, NaN);
    % salinity
    SOCAT_grid.(region{n}).sss_count_nobs = accumarray(subs, SOCAT.(region{n}).salinity, sz, @numel, NaN);
    SOCAT_grid.(region{n}).sss_ave_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @nanmean, NaN);
    SOCAT_grid.(region{n}).sss_std_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @nanstd, NaN);
    SOCAT_grid.(region{n}).sss_max_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @max, NaN);
    SOCAT_grid.(region{n}).sss_min_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @min, NaN);
    
    clear subs sz Xnum Ynum Znum
    
    %% Determine cruise-weighted means and standard deviations
    % If more than one cruise is represented in a given grid cell, replace
    % the unweighted value with a cruise-weighted value
    
    % Pre-allocate with unweighted values
    SOCAT_grid.(region{n}).fco2_ave_wtd  = SOCAT_grid.(region{n}).fco2_ave_unwtd;
    SOCAT_grid.(region{n}).fco2_std_wtd  = SOCAT_grid.(region{n}).fco2_std_unwtd;
    SOCAT_grid.(region{n}).sst_ave_wtd  = SOCAT_grid.(region{n}).sst_ave_unwtd;
    SOCAT_grid.(region{n}).sst_std_wtd  = SOCAT_grid.(region{n}).sst_std_unwtd;
    SOCAT_grid.(region{n}).sss_ave_wtd  = SOCAT_grid.(region{n}).sss_ave_unwtd;
    SOCAT_grid.(region{n}).sss_std_wtd  = SOCAT_grid.(region{n}).sss_std_unwtd;
    SOCAT_grid.(region{n}).fco2_grid_uncert  = nan(size(SOCAT_grid.(region{n}).fco2_std_unwtd));
    % Determine monthly means for individual years
    for a = 1:SOCAT_grid.(region{n}).dim.x
        for b = 1:SOCAT_grid.(region{n}).dim.y
            for c = 1:SOCAT_grid.(region{n}).dim.z
                if SOCAT_grid.(region{n}).count_ncruise(a,b,c) > 1
                    % index to specific grid cell
                    idx = SOCAT.(region{n}).longitude >= SOCAT_grid.(region{n}).lon(a) - 0.125 & ...
                          SOCAT.(region{n}).longitude < SOCAT_grid.(region{n}).lon(a) + 0.125 & ...
                          SOCAT.(region{n}).latitude  >= SOCAT_grid.(region{n}).lat(b) - 0.125 & ...
                          SOCAT.(region{n}).latitude  < SOCAT_grid.(region{n}).lat(b) + 0.125 & ...
                          SOCAT.(region{n}).month  >= SOCAT_grid.(region{n}).month(c) & ...
                          SOCAT.(region{n}).month  < SOCAT_grid.(region{n}).month(c) + 1;
                    cruises = unique(SOCAT.(region{n}).expocode(idx));
                    cruiselist = SOCAT.(region{n}).expocode(idx);
                    fco2 = SOCAT.(region{n}).fCO2(idx);
                    temperature = SOCAT.(region{n}).temperature(idx);
                    salinity = SOCAT.(region{n}).salinity(idx);
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
                    SOCAT_grid.(region{n}).fco2_ave_wtd(a,b,c)  = mean(fco22,'omitnan');
                    SOCAT_grid.(region{n}).fco2_std_wtd(a,b,c)  = mean(std22,'omitnan');
                    if any(startsWith(cruises,'3164'))
                        % checks for grid cells that include moorings; these
                        % are excluded from uncertainty estimates due to
                        % thehigh standard deviation among mooring observations
                        % keyboard
                    else
                        SOCAT_grid.(region{n}).fco2_grid_uncert(a,b,c)  = ...
                            SOCAT_grid.(region{n}).fco2_std_unwtd(a,b,c);
                    end
                    SOCAT_grid.(region{n}).sst_ave_wtd(a,b,c)  = mean(temperature2,'omitnan');
                    SOCAT_grid.(region{n}).sst_std_wtd(a,b,c)  = std(temperature2,[],'omitnan');
                    SOCAT_grid.(region{n}).salinity_ave_wtd(a,b,c)  = mean(salinity2,'omitnan');
                    SOCAT_grid.(region{n}).salinity_std_wtd(a,b,c)  = std(salinity2,[],'omitnan');
                end
            end
        end
    end
    
    clear a b c idx cruises cruiselist fco2 temperature salinity
    clear fco22 std22 temperature2 salinity2 k cruiseidx
    
    %% Count number of months with at least one observation in each grid cell
    SOCAT_grid.(region{n}).num_months = ...
        sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_unwtd),3);
    
    %% Determine area of each grid cell
    SOCAT_grid.(region{n}).area_km2 = ...
    (((repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1) + 0.125) - ...
        (repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1) - 0.125)) .* 110.574) .* ... % latitude distance
    (((repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y) + 0.125) - ...
        (repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y) - 0.125)) .* ...
        111.320.*cosd(repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1))); % longitude distance
    
    %% Determine sea fraction of each grid cell
    load('Data/ETOPO2.mat');
    % limit to LME in question
    lonidx = ETOPO2.lon >= SOCAT_grid.(region{n}).lim.lonmin - 360 & ...
        ETOPO2.lon <= SOCAT_grid.(region{n}).lim.lonmax - 360;
    latidx = ETOPO2.lat >= SOCAT_grid.(region{n}).lim.latmin & ...
        ETOPO2.lat <= SOCAT_grid.(region{n}).lim.latmax;
    ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx);
    ETOPO2.lon = ETOPO2.lon(lonidx);
    ETOPO2.lat = ETOPO2.lat(latidx);
    % define points as land (0) or sea (1)
    ETOPO2.sea = ETOPO2.bottomdepth > 0;
    % determine percentage sea in each RFR-CCS grid cell
    SOCAT_grid.(region{n}).percent_sea = nan(size(SOCAT_grid.(region{n}).lat));
    for a = 1:length(SOCAT_grid.(region{n}).lon)
        for b = 1:length(SOCAT_grid.(region{n}).lat)
            lonidx = find(ETOPO2.lon >= (SOCAT_grid.(region{n}).lon(a)-360)-0.125 & ...
                      ETOPO2.lon < (SOCAT_grid.(region{n}).lon(a)-360)+0.125);
            latidx = find(ETOPO2.lat >= SOCAT_grid.(region{n}).lat(b)-0.125 & ...
                      ETOPO2.lat < SOCAT_grid.(region{n}).lat(b)+0.125);
            SOCAT_grid.(region{n}).percent_sea(a,b) = ...
                sum(sum(ETOPO2.sea(lonidx,latidx)))./...
                (size(ETOPO2.sea(lonidx,latidx),1)*size(ETOPO2.sea(lonidx,latidx),2));
        end
    end
    clear ETOPO2 lonidx latidx a b
    
    %% Plot the percentage of grid cells with data
    figure; worldmap([SOCAT_grid.(region{n}).lim.latmin ...
                      SOCAT_grid.(region{n}).lim.latmax],...
                     [SOCAT_grid.(region{n}).lim.lonmin ...
                      SOCAT_grid.(region{n}).lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1),...
            repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y),...
            SOCAT_grid.(region{n}).num_months);
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

    % save figure
    exportgraphics(gcf,['Figures/' region{n} '_obs.png']);
    close

    % clean up
    clear land c mycolormap

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
    [yf,yr,x] = leastsq(SOCAT_grid.(region{n}).month,...
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

    %% Plot detrended gridded mean pCO2
    figure; worldmap([SOCAT_grid.(region{n}).lim.latmin ...
                      SOCAT_grid.(region{n}).lim.latmax],...
                     [SOCAT_grid.(region{n}).lim.lonmin ...
                      SOCAT_grid.(region{n}).lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1),...
            repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y),...
            mean(SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,3,'omitnan'));
    land = shaperead('landareas', 'UseGeoCoords',true);
    geoshow(land,'FaceColor',rgb('grey'));
    c=colorbar;
    colormap(parula(20));
    caxis([295 475]);
    c.TickLength = 0;
    c.Label.String = 'Surface {\itf}CO_{2} (\muatm)';
    cbarrow;

    % save figure
    exportgraphics(gcf,['Figures/' region{n} '_pCO2.png']);
    close

    % clean up
    clear land c mycolormap

end

% Save gridded pco2 data for all LMEs
save('gridded_pco2','SOCAT_grid','-v7.3');

clear n