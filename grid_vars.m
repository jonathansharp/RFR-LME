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
    lim.latmin = round(min(SOCAT.(region{n}).latitude),0);
    lim.latmax = round(max(SOCAT.(region{n}).latitude),0);
    lim.lonmin = round(min(SOCAT.(region{n}).longitude),0);
    lim.lonmax = round(max(SOCAT.(region{n}).longitude),0);
    lim.monthmin = 1;
    lim.monthmax = 288;
    
    %% Create 0.25 x 0.25 degree grid
%     [SOCAT.([region{n} '_grid']).lat,SOCAT.([region{n} '_grid']).lon,SOCAT.([region{n} '_grid']).month] = ...
%         meshgrid(lim.latmin+0.125:0.25:lim.latmax,lim.lonmin+0.125:0.25:lim.lonmax,lim.monthmin-0.5:1:lim.monthmax-0.5);
%     dim.x = size(SOCAT.([region{n} '_grid']).lon,1);
%     dim.y = size(SOCAT.([region{n} '_grid']).lat,2);
%     dim.z = size(SOCAT.([region{n} '_grid']).month,3);

    SOCAT.([region{n} '_grid']).lon = [lim.lonmin+0.125:0.25:lim.lonmax]';
    dim.x = length(SOCAT.([region{n} '_grid']).lon);
    SOCAT.([region{n} '_grid']).lat = [lim.latmin+0.125:0.25:lim.latmax]';
    dim.y = length(SOCAT.([region{n} '_grid']).lat);
    SOCAT.([region{n} '_grid']).month = [lim.monthmin-0.5:1:lim.monthmax-0.5]';
    dim.z = length(SOCAT.([region{n} '_grid']).month);

    %% Add time variables
    % year
    SOCAT.([region{n} '_grid']).year = repelem(1998:2021,12)';
    % month of year
    SOCAT.([region{n} '_grid']).month_of_year = repmat(1:12,1,24)';
    % date
%     SOCAT.([region{n} '_grid']).date = ...
%         repmat(permute(...
%         datenum([squeeze(SOCAT.([region{n} '_grid']).year(1,1,:)) ...
%         squeeze(SOCAT.([region{n} '_grid']).month(1,1,:)+0.5) ...
%         repmat(15,dim.z,1)]),[3 2 1]),dim.x,dim.y,1);

    %% Determine bin number of each data point
    [~,~,Xnum] = histcounts(SOCAT.(region{n}).longitude,lim.lonmin:0.25:lim.lonmax);
    [~,~,Ynum] = histcounts(SOCAT.(region{n}).latitude,lim.latmin:0.25:lim.latmax);
    [~,~,Znum] = histcounts(SOCAT.(region{n}).month_since_1998,lim.monthmin-1:lim.monthmax);
    
    %% Accumulate 3D grid by applying function to SOCAT values with bin numbers that match grid cells
    subs = [Xnum, Ynum, Znum];
    sz = [dim.x,dim.y,dim.z];
    % cruises
    SOCAT.([region{n} '_grid']).count_ncruise = accumarray(subs, SOCAT.(region{n}).cruise, sz, @(x) numel(unique(x)), NaN);
    % fCO2
    SOCAT.([region{n} '_grid']).fco2_count_nobs = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @numel, NaN);
    SOCAT.([region{n} '_grid']).fco2_ave_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @nanmean, NaN);
    SOCAT.([region{n} '_grid']).fco2_std_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @nanstd, NaN);
    SOCAT.([region{n} '_grid']).fco2_max_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @max, NaN);
    SOCAT.([region{n} '_grid']).fco2_min_unwtd = accumarray(subs, SOCAT.(region{n}).fCO2, sz, @min, NaN);
    % temperature
    SOCAT.([region{n} '_grid']).sst_count_nobs = accumarray(subs, SOCAT.(region{n}).temperature, sz, @numel, NaN);
    SOCAT.([region{n} '_grid']).sst_ave_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @nanmean, NaN);
    SOCAT.([region{n} '_grid']).sst_std_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @nanstd, NaN);
    SOCAT.([region{n} '_grid']).sst_max_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @max, NaN);
    SOCAT.([region{n} '_grid']).sst_min_unwtd = accumarray(subs, SOCAT.(region{n}).temperature, sz, @min, NaN);
    % salinity
    SOCAT.([region{n} '_grid']).sss_count_nobs = accumarray(subs, SOCAT.(region{n}).salinity, sz, @numel, NaN);
    SOCAT.([region{n} '_grid']).sss_ave_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @nanmean, NaN);
    SOCAT.([region{n} '_grid']).sss_std_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @nanstd, NaN);
    SOCAT.([region{n} '_grid']).sss_max_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @max, NaN);
    SOCAT.([region{n} '_grid']).sss_min_unwtd = accumarray(subs, SOCAT.(region{n}).salinity, sz, @min, NaN);
    
    clear subs sz Xnum Ynum Znum
    
    %% Determine cruise-weighted means and standard deviations
    % If more than one cruise is represented in a given grid cell, replace
    % the unweighted value with a cruise-weighted value
    
    % Pre-allocate with unweighted values
    SOCAT.([region{n} '_grid']).fco2_ave_wtd  = SOCAT.([region{n} '_grid']).fco2_ave_unwtd;
    SOCAT.([region{n} '_grid']).fco2_std_wtd  = SOCAT.([region{n} '_grid']).fco2_std_unwtd;
    SOCAT.([region{n} '_grid']).sst_ave_wtd  = SOCAT.([region{n} '_grid']).sst_ave_unwtd;
    SOCAT.([region{n} '_grid']).sst_std_wtd  = SOCAT.([region{n} '_grid']).sst_std_unwtd;
    SOCAT.([region{n} '_grid']).sss_ave_wtd  = SOCAT.([region{n} '_grid']).sss_ave_unwtd;
    SOCAT.([region{n} '_grid']).sss_std_wtd  = SOCAT.([region{n} '_grid']).sss_std_unwtd;
    SOCAT.([region{n} '_grid']).fco2_grid_uncert  = nan(size(SOCAT.([region{n} '_grid']).fco2_std_unwtd));
    % Determine monthly means for individual years
    for a = 1:dim.x
        for b = 1:dim.y
            for c = 1:dim.z
                if SOCAT.([region{n} '_grid']).count_ncruise(a,b,c) > 1
                    % index to specific grid cell
                    idx = SOCAT.(region{n}).longitude >= SOCAT.([region{n} '_grid']).lon(a) - 0.125 & ...
                          SOCAT.(region{n}).longitude < SOCAT.([region{n} '_grid']).lon(a) + 0.125 & ...
                          SOCAT.(region{n}).latitude  >= SOCAT.([region{n} '_grid']).lat(b) - 0.125 & ...
                          SOCAT.(region{n}).latitude  < SOCAT.([region{n} '_grid']).lat(b) + 0.125 & ...
                          SOCAT.(region{n}).month  >= SOCAT.([region{n} '_grid']).month(c) & ...
                          SOCAT.(region{n}).month  < SOCAT.([region{n} '_grid']).month(c) + 1;
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
                    SOCAT.([region{n} '_grid']).fco2_ave_wtd(a,b,c)  = mean(fco22,'omitnan');
                    SOCAT.([region{n} '_grid']).fco2_std_wtd(a,b,c)  = mean(std22,'omitnan');
                    if any(startsWith(cruises,'3164'))
                        % checks for grid cells that include moorings; these
                        % are excluded from uncertainty estimates due to
                        % thehigh standard deviation among mooring observations
                        % keyboard
                    else
                        SOCAT.([region{n} '_grid']).fco2_grid_uncert(a,b,c)  = ...
                            SOCAT.([region{n} '_grid']).fco2_std_unwtd(a,b,c);
                    end
                    SOCAT.([region{n} '_grid']).sst_ave_wtd(a,b,c)  = mean(temperature2,'omitnan');
                    SOCAT.([region{n} '_grid']).sst_std_wtd(a,b,c)  = std(temperature2,[],'omitnan');
                    SOCAT.([region{n} '_grid']).salinity_ave_wtd(a,b,c)  = mean(salinity2,'omitnan');
                    SOCAT.([region{n} '_grid']).salinity_std_wtd(a,b,c)  = std(salinity2,[],'omitnan');
                end
            end
        end
    end
    
    clear a b c idx cruises cruiselist fco2 temperature salinity
    clear fco22 std22 temperature2 salinity2 k cruiseidx
    
    %% Count number of months with at least one observation in each grid cell
    SOCAT.([region{n} '_grid']).num_months = ...
        sum(~isnan(SOCAT.([region{n} '_grid']).fco2_ave_unwtd),3);
    
    %% Determine area of each grid cell
    SOCAT.([region{n} '_grid']).area_km2 = ...
    (((repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1) + 0.125) - ...
        (repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1) - 0.125)) .* 110.574) .* ... % latitude distance
    (((repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y) + 0.125) - ...
        (repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y) - 0.125)) .* ...
        111.320.*cosd(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1))); % longitude distance
    
    %% Determine sea fraction of each grid cell
    load('Data/ETOPO2.mat');
    % limit to LME in question
    lonidx = ETOPO2.lon >= lim.lonmin - 360 & ETOPO2.lon <= lim.lonmax - 360;
    latidx = ETOPO2.lat >= lim.latmin & ETOPO2.lat <= lim.latmax;
    ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx);
    ETOPO2.lon = ETOPO2.lon(lonidx);
    ETOPO2.lat = ETOPO2.lat(latidx);
    % define points as land (0) or sea (1)
    ETOPO2.sea = ETOPO2.bottomdepth > 0;
    % determine percentage sea in each RFR-CCS grid cell
    SOCAT.([region{n} '_grid']).percent_sea = nan(size(SOCAT.([region{n} '_grid']).lat));
    for a = 1:length(SOCAT.([region{n} '_grid']).lon)
        for b = 1:length(SOCAT.([region{n} '_grid']).lat)
            lonidx = find(ETOPO2.lon >= (SOCAT.([region{n} '_grid']).lon(a)-360)-0.125 & ...
                      ETOPO2.lon < (SOCAT.([region{n} '_grid']).lon(a)-360)+0.125);
            latidx = find(ETOPO2.lat >= SOCAT.([region{n} '_grid']).lat(b)-0.125 & ...
                      ETOPO2.lat < SOCAT.([region{n} '_grid']).lat(b)+0.125);
            SOCAT.([region{n} '_grid']).percent_sea(a,b) = ...
                sum(sum(ETOPO2.sea(lonidx,latidx)))./...
                (size(ETOPO2.sea(lonidx,latidx),1)*size(ETOPO2.sea(lonidx,latidx),2));
        end
    end
    clear ETOPO2 lonidx latidx a b
    
    %% Plot the percentage of grid cells with data
    figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
            repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
            SOCAT.([region{n} '_grid']).num_months);
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
    area_weights = SOCAT.([region{n} '_grid']).area_km2;
    area_weights = repmat(area_weights,1,1,size(SOCAT.([region{n} '_grid']).fco2_ave_wtd,3));
    area_weights(isnan(SOCAT.([region{n} '_grid']).fco2_ave_wtd)) = NaN;
    % Calculate area-weighted domain mean
    SOCAT.([region{n} '_grid']).fco2_dom_mean = ...
        squeeze(sum(sum(SOCAT.([region{n} '_grid']).fco2_ave_wtd.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    SOCAT.([region{n} '_grid']).fco2_dom_mean(SOCAT.([region{n} '_grid']).fco2_dom_mean == 0) = NaN;
    % Fit trend to area weighted domain mean
    [yf,yr,x] = leastsq(SOCAT.([region{n} '_grid']).month,...
        SOCAT.([region{n} '_grid']).fco2_dom_mean,0,0,0);
    % Remove difference from mean for each month
    for m = 1:length(SOCAT.([region{n} '_grid']).month)
        SOCAT.([region{n} '_grid']).fco2_ave_wtd_detrend(:,:,m) = ...
            SOCAT.([region{n} '_grid']).fco2_ave_wtd(:,:,m) + ...
            (mean(yf,'omitnan') - yf(m,:));
    end
    % Calculate area-weighted detrended domain mean
    SOCAT.([region{n} '_grid']).fco2_dom_mean_detrend = ...
        squeeze(sum(sum(SOCAT.([region{n} '_grid']).fco2_ave_wtd_detrend.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    SOCAT.([region{n} '_grid']).fco2_dom_mean_detrend(SOCAT.([region{n} '_grid']).fco2_dom_mean_detrend == 0) = NaN;

    %% Plot detrended gridded mean pCO2
    figure; worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(SOCAT.([region{n} '_grid']).lat',dim.x,1),...
            repmat(SOCAT.([region{n} '_grid']).lon,1,dim.y),...
            mean(SOCAT.([region{n} '_grid']).fco2_ave_wtd_detrend,3,'omitnan'));
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

clear n