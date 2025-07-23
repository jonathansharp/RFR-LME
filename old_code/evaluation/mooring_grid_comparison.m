% This script compares RFR-LME pCO2 grids constructed with and without
% mooring observations

% Load gridded mooring observations
load('Data/socat_gridded_2022_moorings_only','SOCAT_grid');
% define mooring indices by grid cells that include data
[a,b] = find(~isnan(mean(SOCAT_grid.fco2_ave_wtd,3,'omitnan')));

% extract mooring timeseries
cnt=1; % set counter to one
for m = 1:length(a)
    mooring.(['m' num2str(m)]).fCO2 = ...
        squeeze(SOCAT_grid.fco2_ave_wtd(a(m),b(m),:));
    mooring.(['m' num2str(m)]).SST = ...
        squeeze(SOCAT_grid.sst_ave_wtd(a(m),b(m),:));
    mooring.(['m' num2str(m)]).SSS = ...
        squeeze(SOCAT_grid.sss_ave_wtd(a(m),b(m),:));
    C=CO2SYS(2400,mooring.(['m' num2str(m)]).fCO2,1,5,mooring.(['m' num2str(m)]).SSS,...
        mooring.(['m' num2str(m)]).SST,NaN,0,NaN,0,0,0,0,1,10,1,2,2);
    pCO2 = C(:,4); pCO2(pCO2==-999)=NaN;
    mooring.(['m' num2str(m)]).pCO2 = pCO2;
    mooring.(['m' num2str(m)]).lon = SOCAT_grid.lon(a(m));
    mooring.(['m' num2str(m)]).lat = SOCAT_grid.lat(b(m));
    mooring.(['m' num2str(m)]).time = datenum(1998,SOCAT_grid.month+0.5,15);
    % average together the two WHOTS grid cells (since the lat/lon changes a bit
    if mooring.(['m' num2str(m)]).lon == 202.125 && ...
            mooring.(['m' num2str(m)]).lat > 22.5 && ...
            mooring.(['m' num2str(m)]).lat < 23
        temp(:,cnt) = mooring.(['m' num2str(m)]).pCO2;
        WHOTSm(cnt) = m;
        cnt=cnt+1;
    end
    % remove moorings that aren't in LMEs or don't have enough obs
%     if sum(~isnan(mooring.(['m' num2str(m)]).pCO2)) < 36
%         mooring = rmfield(mooring,(['m' num2str(m)]));
%     end
end
mooring.(['m' num2str(WHOTSm(1))]).pCO2 = mean(temp,2,'omitnan');
mooring = rmfield(mooring,(['m' num2str(WHOTSm(2))]));
moor_nums = fieldnames(mooring);
clear a b m cnt temp WHOTSm % clean up

% this script defines the bounds of the eleven LMEs
define_regions_eiwg

%% loop through each region
for n = 7%1:length(region)

    % load datasets
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    w_moor = load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    no_moor = load(['Data/' region{n} '/ML_fCO2_no_moorings'],'OAI_grid');

    % calculates differences and their variability
    diffs = w_moor.OAI_grid.(region{n}).pCO2 - no_moor.OAI_grid.(region{n}).pCO2;
    diffs_rmsd = sqrt(mean(diffs.^2,3,'omitnan'));
    diffs_abs_avg = mean(abs(diffs),3,'omitnan');

    % plot average absolute differences
    figure('visible','on');
    set(gca,'fontsize',14);
    worldmap([w_moor.OAI_grid.(region{n}).lim.latmin ...
        w_moor.OAI_grid.(region{n}).lim.latmax],...
        [w_moor.OAI_grid.(region{n}).lim.lonmin ...
        w_moor.OAI_grid.(region{n}).lim.lonmax]);
    contourfm(w_moor.OAI_grid.(region{n}).lat,w_moor.OAI_grid.(region{n}).lon,...
        diffs_abs_avg',0:0.5:20,'LineStyle','none');
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    plot_land('map');
    c=colorbar;
    colormap(cmocean('amp'));
    caxis([0 20]);
    c.TickLength = 0;
    c.Label.String = 'Mean |\Delta{\itp}CO_{2}| (\muatm)';
    cbarrow('up');
    mlabel off;
    plabel off;
    % save figure (without labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_abs_avg.png'],'-transparent');
    % plot pCO2 from moorings
    for m = 1:length(moor_nums)
        if inpolygon(mooring.(moor_nums{m}).lon,mooring.(moor_nums{m}).lat,...
                convert_lon(lme_shape(lme_idx.(region{n})).X),lme_shape(lme_idx.(region{n})).Y)
            scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
                100,'pentagram','MarkerEdgeColor','k','MarkerFaceColor','y');
        end
    end
    % save figure (with labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_abs_avg_w_stars.png'],'-transparent');
    close
    
    % plot variability of differences
    figure('visible','on');
    set(gca,'fontsize',14);
    worldmap([w_moor.OAI_grid.(region{n}).lim.latmin ...
        w_moor.OAI_grid.(region{n}).lim.latmax],...
        [w_moor.OAI_grid.(region{n}).lim.lonmin ...
        w_moor.OAI_grid.(region{n}).lim.lonmax]);
    contourfm(w_moor.OAI_grid.(region{n}).lat,w_moor.OAI_grid.(region{n}).lon,...
        diffs_rmsd',0:0.5:20,'LineStyle','none');
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    plot_land('map');
    c=colorbar;
    colormap(cmocean('amp'));
    caxis([0 20]);
    c.TickLength = 0;
    c.Label.String = '\Delta{\itp}CO_{2} RMSD (\muatm)';
    cbarrow('up');
    mlabel off;
    plabel off;
    % save figure (without labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_diffs.png'],'-transparent');
    % plot pCO2 from moorings
    for m = 1:length(moor_nums)
        if inpolygon(mooring.(moor_nums{m}).lon,mooring.(moor_nums{m}).lat,...
                convert_lon(lme_shape(lme_idx.(region{n})).X),lme_shape(lme_idx.(region{n})).Y)
            scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
                100,'pentagram','MarkerEdgeColor','k','MarkerFaceColor','y');
        end
    end
    % save figure (with labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_diffs_w_stars.png'],'-transparent');
    close

    % calculate amplitude and differences
    w_moor_pco2_clim = nan(w_moor.OAI_grid.(region{n}).dim.x,w_moor.OAI_grid.(region{n}).dim.y,12);
    no_moor_pco2_clim = nan(no_moor.OAI_grid.(region{n}).dim.x,no_moor.OAI_grid.(region{n}).dim.y,12);
    for m=1:12
        w_moor_pco2_clim(:,:,m) = mean(w_moor.OAI_grid.(region{n}).pCO2(:,:,m:12:end),3,'omitnan');
        no_moor_pco2_clim(:,:,m) = mean(no_moor.OAI_grid.(region{n}).pCO2(:,:,m:12:end),3,'omitnan');
    end
    w_moor_amp = max(w_moor_pco2_clim,[],3) - min(w_moor_pco2_clim,[],3);
    no_moor_amp = max(no_moor_pco2_clim,[],3) - min(no_moor_pco2_clim,[],3);
    amp_diff = w_moor_amp - no_moor_amp;

    % plot variability of differences
    figure('visible','on');
    set(gca,'fontsize',14);
    worldmap([w_moor.OAI_grid.(region{n}).lim.latmin ...
        w_moor.OAI_grid.(region{n}).lim.latmax],...
        [w_moor.OAI_grid.(region{n}).lim.lonmin ...
        w_moor.OAI_grid.(region{n}).lim.lonmax]);
    contourfm(w_moor.OAI_grid.(region{n}).lat,w_moor.OAI_grid.(region{n}).lon,...
        amp_diff',-20:1:20,'LineStyle','none');
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    caxis([-10 10]);
    plot_land('map');
    c=colorbar;
    colormap(cmocean('balance','pivot',0));
    c.TickLength = 0;
    c.Label.String = '\Delta{\itp}CO_{2} Amp. (\muatm)';
    cbarrow;
    mlabel off;
    plabel off;
    % save figure (without labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_amp_diffs.png'],'-transparent');
    % plot pCO2 from moorings
    for m = 1:length(moor_nums)
        if inpolygon(mooring.(moor_nums{m}).lon,mooring.(moor_nums{m}).lat,...
            convert_lon(lme_shape(lme_idx.(region{n})).X),lme_shape(lme_idx.(region{n})).Y)
            scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
                100,'pentagram','MarkerEdgeColor','k','MarkerFaceColor','y');
        end
    end
    % save figure (with labels)
    if ~isfolder('Figures'); mkdir('Figures'); end
    export_fig(gcf,['Figures/' region{n} '_pCO2_moor_amp_diffs_w_stars.png'],'-transparent');
    close

end


%% zoomed in figure

% NE region
n=7;

% load datasets
load(['Data/' (region{n}) '/gridded_pco2'],'SOCAT_grid');
w_moor = load(['Data/' (region{n}) '/ML_fCO2'],'OAI_grid');
no_moor = load(['Data/' (region{n}) '/ML_fCO2_no_moorings'],'OAI_grid');

% calculates differences and their variability
diffs = w_moor.OAI_grid.(region{n}).pCO2 - no_moor.OAI_grid.(region{n}).pCO2;
diffs_rmsd = sqrt(mean(diffs.^2,3,'omitnan'));
diffs_abs_avg = mean(abs(diffs),3,'omitnan');
num=10; highest = maxk(diffs_rmsd(:),num);
lon_idx = nan(num,1); lat_idx = nan(num,1);
for nn = 1:num
    [lon_idx(nn),lat_idx(nn)] = find(diffs_rmsd==highest(nn));
end

% plot variability of differences
figure('visible','on');
set(gca,'fontsize',14);
worldmap([41 44],[288 291]);
% contourfm(w_moor.OAI_grid.(region{n}).lat,w_moor.OAI_grid.(region{n}).lon,...
%     diffs_rmsd',0:0.5:20,'LineStyle','none');
pcolorm(repmat(w_moor.OAI_grid.(region{n}).lat-.125,1,length(w_moor.OAI_grid.(region{n}).lon)),...
    repmat(w_moor.OAI_grid.(region{n}).lon'-.125,length(w_moor.OAI_grid.(region{n}).lat),1),...
    diffs_rmsd');
scatterm(w_moor.OAI_grid.(region{n}).lat(lat_idx),w_moor.OAI_grid.(region{n}).lon(lon_idx),...
    300,'k.');
% plot land
plot_land('map');
c=colorbar;
clim([0 20]);
colormap(cmocean('amp'));
c.TickLength = 0;
c.Label.String = '\Delta{\itp}CO_{2} RMSD (\muatm)';
cbarrow('up');
mlabel off;
plabel off;
% save figure (without labels)
if ~isfolder('Figures'); mkdir('Figures'); end
export_fig(gcf,['Figures/' (region{n}) '_pCO2_moor_diffs_zoomed.png'],'-transparent');
% plot pCO2 from moorings
for m = 1:length(moor_nums)
    if inpolygon(mooring.(moor_nums{m}).lon,mooring.(moor_nums{m}).lat,...
            convert_lon(lme_shape(lme_idx.(region{n})).X),lme_shape(lme_idx.(region{n})).Y)
        scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
            100,'pentagram','MarkerEdgeColor','k','MarkerFaceColor','y');
    end
end
% save figure (with labels)
if ~isfolder('Figures'); mkdir('Figures'); end
export_fig(gcf,['Figures/' (region{n}) '_pCO2_moor_diffs_w_stars_zoomed.png'],'-transparent');
close

% calculate amplitude and differences
w_moor_pco2_clim = nan(w_moor.OAI_grid.(region{n}).dim.x,w_moor.OAI_grid.(region{n}).dim.y,12);
no_moor_pco2_clim = nan(no_moor.OAI_grid.(region{n}).dim.x,no_moor.OAI_grid.(region{n}).dim.y,12);
for m=1:12
    w_moor_pco2_clim(:,:,m) = mean(w_moor.OAI_grid.(region{n}).pCO2(:,:,m:12:end),3,'omitnan');
    no_moor_pco2_clim(:,:,m) = mean(no_moor.OAI_grid.(region{n}).pCO2(:,:,m:12:end),3,'omitnan');
end
w_moor_amp = max(w_moor_pco2_clim,[],3) - min(w_moor_pco2_clim,[],3);
no_moor_amp = max(no_moor_pco2_clim,[],3) - min(no_moor_pco2_clim,[],3);
amp_diff = w_moor_amp - no_moor_amp;
num=10; highest = maxk(abs(amp_diff(:)),num);
lon_idx = nan(num,1); lat_idx = nan(num,1);
for nn = 1:num
    [lon_idx(nn),lat_idx(nn)] = find(abs(amp_diff)==highest(nn));
end

% plot amplitude differences
figure('visible','on');
set(gca,'fontsize',14);
worldmap([41 44],[288 291]);
% contourfm(w_moor.OAI_grid.(region{n}).lat,w_moor.OAI_grid.(region{n}).lon,...
%     diffs_rmsd',0:0.5:20,'LineStyle','none');
pcolorm(repmat(w_moor.OAI_grid.(region{n}).lat-.125,1,length(w_moor.OAI_grid.(region{n}).lon)),...
    repmat(w_moor.OAI_grid.(region{n}).lon'-.125,length(w_moor.OAI_grid.(region{n}).lat),1),...
    amp_diff');
scatterm(w_moor.OAI_grid.(region{n}).lat(lat_idx),w_moor.OAI_grid.(region{n}).lon(lon_idx),...
    300,'k.');
% plot land
plot_land('map');
c=colorbar;
clim([-10 10]);
colormap(cmocean('balance','pivot',0));
c.TickLength = 0;
c.Label.String = '\Delta{\itp}CO_{2} Amp. (\muatm)';
mlabel off;
plabel off;
% save figure (without labels)
if ~isfolder('Figures'); mkdir('Figures'); end
export_fig(gcf,['Figures/' (region{n}) '_pCO2_moor_diffs_amp_zoomed.png'],'-transparent');
% plot pCO2 from moorings
for m = 1:length(moor_nums)
    if inpolygon(mooring.(moor_nums{m}).lon,mooring.(moor_nums{m}).lat,...
            convert_lon(lme_shape(lme_idx.(region{n})).X),lme_shape(lme_idx.(region{n})).Y)
        scatterm(mooring.(moor_nums{m}).lat,mooring.(moor_nums{m}).lon,...
            100,'pentagram','MarkerEdgeColor','k','MarkerFaceColor','y');
    end
end
% save figure (with labels)
if ~isfolder('Figures'); mkdir('Figures'); end
export_fig(gcf,['Figures/' (region{n}) '_pCO2_moor_diffs_amp_w_stars_zoomed.png'],'-transparent');
close
