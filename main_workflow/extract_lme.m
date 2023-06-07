% Extract regional SOCAT data
% 
% This script grids observations of fCO2 and ancillary variables from five
% US large Marine Ecosystems (Alaska, California Current, Insular Pacific /
% Hawaii, Gulf of Mexico / Caribbean, and US East Coast) assembled via
% extractions from SOCATv2022 defined by latitude and longitude bounds.
% 
% Written by J.D. Sharp: 7/26/22
% Last updated by J.D. Sharp: 5/15/22
% 

%% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

%% extract each LME from large grid
for n = 1:length(region)

    %% load SOCAT grid
    load('Data/socat_gridded','SOCAT_grid');

    %% display status
    disp(['Extracting ' region{n} ' LME from SOCAT grid']);

    %% remove observations outside general LME limits
    % determine geographic indices
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X)';
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
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

    %% visualize number of observations
%     time = datetime(SOCAT.year,...
%         SOCAT.month,SOCAT.day);
%     figure('visible','off');
%     histogram(time);
%     ylabel('Number of observations');
%     xlabel('Year');
%     clear time
% 
%     % save figure
%     if ~isfolder(['Figures/' region{n}]); mkdir(['Figures/' region{n}]); end
%     exportgraphics(gcf,['Figures/' region{n} '/hist.png']);
%     close

    plot_regional_gif(SOCAT_grid.(region{n}).lim,...
        SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,...
        SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,parula(20),'fCO2_obs',...
        'Surface {\itf}CO_{2} (\muatm)',SOCAT_grid.(region{n}).year,...
        SOCAT_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));

    %% Plot the percentage of grid cells with data
    figure('visible','off');
    worldmap([SOCAT_grid.(region{n}).lim.latmin ...
        SOCAT_grid.(region{n}).lim.latmax],...
       [SOCAT_grid.(region{n}).lim.lonmin ...
        SOCAT_grid.(region{n}).lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(SOCAT_grid.(region{n}).lat',SOCAT_grid.(region{n}).dim.x,1),...
            repmat(SOCAT_grid.(region{n}).lon,1,SOCAT_grid.(region{n}).dim.y),...
            SOCAT_grid.(region{n}).num_months);
    plot_land('map');
    c=colorbar;
    mycolormap = jet(21);
    mycolormap(1,:) = 1;
    colormap(mycolormap);
    caxis([-1 41]);
    c.TickLength = 0;
    c.Label.String = 'Number of Months Represented';
    cbarrow('up');

    % save figure
    if ~isfolder(['Figures/' region{n}]); mkdir(['Figures/' region{n}]); end
    exportgraphics(gcf,['Figures/' region{n} '_num_obs.png']);
    close

    % clean up
    clear c mycolormap

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

    %% Plot detrended gridded mean pCO2
    plot_temporal_mean(SOCAT_grid.(region{n}).lim,...
        SOCAT_grid.(region{n}).dim,SOCAT_grid.(region{n}).lat,...
        SOCAT_grid.(region{n}).lon,SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,...
        parula(20),'fCO2_obs','Surface {\itf}CO_{2} (\muatm)',region{n});

    % Save gridded pco2 data
    if ~isfolder(['Data/' region{n}]); mkdir(['Data/' region{n}]); end
    save(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid','-v7.3');

    % clean up
    clear SOCAT_grid

end

%% Plot number of all grid cells with data
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
mycolormap = jet(21);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-1 41]);
c.TickLength = 0;
c.Label.String = 'Total Months Represented';
cbarrow('right');
% plot background
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,SOCAT_grid.num_months');
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,...
        SOCAT_grid.(region{n}).num_months');
    clear SOCAT_grid
end
% plot borders around regions
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% plot land
plot_land('map');
mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/num_obs.png');
close
% clean up
clear n h r tmp_lon c mycolormap 

%% Plot number of climatological grid cells with data
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
mycolormap = jet(13);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-0.5 12.5]);
c.TickLength = 0;
c.Label.String = 'Months of Year Represented';
% plot background
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,SOCAT_grid.num_months_clim');
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,...
        SOCAT_grid.(region{n}).num_months_clim');
    clear SOCAT_grid
end
% plot borders around regions
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% plot land
plot_land('map');
mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/num_obs_clim.png');
close
% clean up
clear n h r tmp_lon c mycolormap 

%% Plot all detrended gridded mean pCO2
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
% title('DJF');
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(parula(20));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Surface {\itf}CO_{2} (\muatm)';
cbarrow;
% plot background
% z = mean(SOCAT_grid.fco2_ave_wtd_detrend,3,'omitnan')';
z = mean(cat(3,SOCAT_grid.fco2_ave_wtd_detrend(:,:,9:12:end),...
               SOCAT_grid.fco2_ave_wtd_detrend(:,:,10:12:end),...
               SOCAT_grid.fco2_ave_wtd_detrend(:,:,11:12:end)),...
               3,'omitnan')';
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,z);
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
%     z = mean(SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,3,'omitnan')';
    z = mean(cat(3,SOCAT_grid.(region{n}).fco2_ave_wtd_detrend(:,:,9:12:end),...
               SOCAT_grid.(region{n}).fco2_ave_wtd_detrend(:,:,10:12:end),...
               SOCAT_grid.(region{n}).fco2_ave_wtd_detrend(:,:,11:12:end)),...
               3,'omitnan')';
    pcolorm(SOCAT_grid.(region{n}).lat,SOCAT_grid.(region{n}).lon,z);
    clear SOCAT_grid
end
% plot borders around regions
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% figure properties
plot_land('map');
mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,'Figures/full/fCO2_obs_SON.png');
close
% clean up
clear n z h r c tmp_lon
