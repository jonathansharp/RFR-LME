% Compare to mooring observations
% 
% This script compares fCO2 predictions in US LMEs to gridded fCO2 at
% mooring observation sites, obtained from SOCAT
% 
% Written by J.D. Sharp: 8/1/23
% Last updated by J.D. Sharp: 8/1/23

%% Load gridded mooring observations
load('Data/socat_gridded_2023_moorings_only','SOCAT_grid');
% define mooring indices
[a,b] = find(~isnan(mean(SOCAT_grid.fco2_ave_wtd,3,'omitnan')));
% extract mooring timeseries
for m = 1:length(a)
    mooring.(['m' num2str(m) '_fCO2']) = ...
        squeeze(SOCAT_grid.fco2_ave_wtd(a(m),b(m),:));
end

%% Pre-allocate LME RFR grid
US_LME_RFR.lim = SOCAT_grid.lim;
US_LME_RFR.dim = SOCAT_grid.dim;
US_LME_RFR.lon = SOCAT_grid.lon;
US_LME_RFR.lat = SOCAT_grid.lat;
US_LME_RFR.month = SOCAT_grid.month;
US_LME_RFR.year = SOCAT_grid.year;
US_LME_RFR.fCO2_no_moorings = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
US_LME_RFR.fCO2_full = nan(US_LME_RFR.dim.x,US_LME_RFR.dim.y,US_LME_RFR.dim.z);
clear SOCAT_grid m

%% Add mooring-excluded fCO2 to full grid
define_regions_eiwg
for n = 1:length(region)
    % load gridded fCO2
    load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');
    % Add fCO2 to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);
    US_LME_RFR.fCO2_no_moorings(idx_lon,idx_lat,:) = OAI_grid.(region{n}).fCO2;
%     US_LME_RFR.ufco2(idx_lon,idx_lat,:) = OAI_grid.(region{n}).ufCO2;
    % clean up
    clear idx_lon idx_lat OAI_grid n
end
%% Add full fCO2 to full grid
define_regions_eiwg
for n = 1:length(region)
    % load gridded fCO2
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % Add fCO2 to full grid
    idx_lon = US_LME_RFR.lon >= min(OAI_grid.(region{n}).lon) & ...
        US_LME_RFR.lon <= max(OAI_grid.(region{n}).lon);
    idx_lat = US_LME_RFR.lat >= min(OAI_grid.(region{n}).lat) & ...
        US_LME_RFR.lat <= max(OAI_grid.(region{n}).lat);
    US_LME_RFR.fCO2_full(idx_lon,idx_lat,:) = OAI_grid.(region{n}).fCO2;
%     US_LME_RFR.ufco2(idx_lon,idx_lat,:) = OAI_grid.(region{n}).ufCO2;
    % clean up
    clear idx_lon idx_lat OAI_grid n
end

%% plot mooring timeseries comparisons
clrs = cbrewer('qual','Dark2',3);
for m = 1:length(a)
    if ~isnan(US_LME_RFR.fCO2_full(a(m),b(m),1)) && ...
            sum(~isnan(mooring.(['m' num2str(m) '_fCO2']))) > 36
        % create plot
        figure; hold on
        set(gcf,'position',[10 10 800 400]);
        title([num2str(US_LME_RFR.lat(b(m))) char(176) 'N, ' ...
            num2str(US_LME_RFR.lon(a(m))) char(176) 'E']);
        p1=plot(US_LME_RFR.month,squeeze(US_LME_RFR.fCO2_no_moorings(a(m),b(m),:)),...
            'color',clrs(1,:),'linewidth',2);
        p2=plot(US_LME_RFR.month,squeeze(US_LME_RFR.fCO2_full(a(m),b(m),:)),...
            ':','color',clrs(2,:),'linewidth',2);
        scatter(US_LME_RFR.month,mooring.(['m' num2str(m) '_fCO2']),...
            'MarkerFaceColor',clrs(3,:),'MarkerEdgeColor',clrs(3,:));
        % add stats
        del_no_moorings = mooring.(['m' num2str(m) '_fCO2']) - ...
            squeeze(US_LME_RFR.fCO2_no_moorings(a(m),b(m),:));
        del_full = mooring.(['m' num2str(m) '_fCO2']) - ...
            squeeze(US_LME_RFR.fCO2_full(a(m),b(m),:));
        avg_del_no_moorings = mean(del_no_moorings,'omitnan');
        avg_del_full = mean(del_full,'omitnan');
        std_del_no_moorings = std(del_no_moorings,[],'omitnan');
        std_del_full = std(del_full,[],'omitnan');
        legend([p1 p2],{['No Moorings: ' num2str(round(avg_del_no_moorings,1)) ' ' char(177) ' ' ...
            num2str(round(std_del_no_moorings,1))],['Moorings: ' ...
            num2str(round(avg_del_full,1)) ' ' char(177) ' ' num2str(round(std_del_full,1))]},...
            'location','northwest','fontsize',12);
        % export figure
        exportgraphics(gcf,['Figures/moorings/m' num2str(m) '.png']);
        close
    end
end

%% plot mooring map
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(parula);
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itf}_{CO2} (\muatm)';
cbarrow;
% plot regions
for n = 1:length(region)
    vars_grid = load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');
    z = mean(vars_grid.OAI_grid.(region{n}).fCO2,3,'omitnan')';
    contourfm(vars_grid.OAI_grid.(region{n}).lat,vars_grid.OAI_grid.(region{n}).lon,...
        z,295:(475-295)/200:475,'LineStyle','none');
    clear vars_grid z
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
% plot mooring locations
for m = 1:length(a)
    if ~isnan(US_LME_RFR.fCO2_full(a(m),b(m),1)) && ...
            sum(~isnan(mooring.(['m' num2str(m) '_fCO2']))) > 36
        scatterm(US_LME_RFR.lat(b(m)),US_LME_RFR.lon(a(m)),80,'.r')
    end
end
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,['Figures/full/fCO2_with_moorings.png']);
close
