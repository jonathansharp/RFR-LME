% Evaluate US LME OA indicators at specific locations

%% load US LME data
date = '20-May-2023';
LME_RFR = netcdfreader(['Data/US_LME_RFR_Inds_' date '.nc']);
LME_RFR_u = netcdfreader(['Data/US_LME_RFR_Inds_Uncer_' date '.nc']);

%% extract time series from RFR data
Time = double(LME_RFR.Time);
% Northern Gulf of Alaska time series
NGoA.lon = 214.125;
NGoA.lat = 59.625;
NGoA.lon_idx = find(abs(LME_RFR.Lon - NGoA.lon) == min(abs(LME_RFR.Lon - NGoA.lon)));
NGoA.lat_idx = find(abs(LME_RFR.Lat - NGoA.lat) == min(abs(LME_RFR.Lat - NGoA.lat)));
NGoA.fCO2 = squeeze(LME_RFR.fCO2(NGoA.lon_idx,NGoA.lat_idx,:));
NGoA.ufCO2 = squeeze(LME_RFR_u.ufCO2(NGoA.lon_idx,NGoA.lat_idx,:));
NGoA.pH = squeeze(LME_RFR.pH(NGoA.lon_idx,NGoA.lat_idx,:));
NGoA.upH = squeeze(LME_RFR_u.upH(NGoA.lon_idx,NGoA.lat_idx,:));
% Southern Gulf of Alaska time series
SGoA.lon = 225.375;
SGoA.lat = 55.125;
SGoA.lon_idx = find(abs(LME_RFR.Lon - SGoA.lon) == min(abs(LME_RFR.Lon - SGoA.lon)));
SGoA.lat_idx = find(abs(LME_RFR.Lat - SGoA.lat) == min(abs(LME_RFR.Lat - SGoA.lat)));
SGoA.fCO2 = squeeze(LME_RFR.fCO2(SGoA.lon_idx,SGoA.lat_idx,:));
SGoA.ufCO2 = squeeze(LME_RFR_u.ufCO2(SGoA.lon_idx,SGoA.lat_idx,:));
SGoA.pH = squeeze(LME_RFR.pH(SGoA.lon_idx,SGoA.lat_idx,:));
SGoA.upH = squeeze(LME_RFR_u.upH(SGoA.lon_idx,SGoA.lat_idx,:));
% Northern California Current time series
NCal.lon = 235.125;
NCal.lat = 45.875;
NCal.lon_idx = find(abs(LME_RFR.Lon - NCal.lon) == min(abs(LME_RFR.Lon - NCal.lon)));
NCal.lat_idx = find(abs(LME_RFR.Lat - NCal.lat) == min(abs(LME_RFR.Lat - NCal.lat)));
NCal.fCO2 = squeeze(LME_RFR.fCO2(NCal.lon_idx,NCal.lat_idx,:));
NCal.ufCO2 = squeeze(LME_RFR_u.ufCO2(NCal.lon_idx,NCal.lat_idx,:));
NCal.pH = squeeze(LME_RFR.pH(NCal.lon_idx,NCal.lat_idx,:));
NCal.upH = squeeze(LME_RFR_u.upH(NCal.lon_idx,NCal.lat_idx,:));
% Southern California Current time series
MCal.lon = 235.375;
MCal.lat = 39.125;
MCal.lon_idx = find(abs(LME_RFR.Lon - MCal.lon) == min(abs(LME_RFR.Lon - MCal.lon)));
MCal.lat_idx = find(abs(LME_RFR.Lat - MCal.lat) == min(abs(LME_RFR.Lat - MCal.lat)));
MCal.fCO2 = squeeze(LME_RFR.fCO2(MCal.lon_idx,MCal.lat_idx,:));
MCal.ufCO2 = squeeze(LME_RFR_u.ufCO2(MCal.lon_idx,MCal.lat_idx,:));
MCal.pH = squeeze(LME_RFR.pH(MCal.lon_idx,MCal.lat_idx,:));
MCal.upH = squeeze(LME_RFR_u.upH(MCal.lon_idx,MCal.lat_idx,:));
% Southern California Current time series
SCal.lon = 241.625;
SCal.lat = 32.875;
SCal.lon_idx = find(abs(LME_RFR.Lon - SCal.lon) == min(abs(LME_RFR.Lon - SCal.lon)));
SCal.lat_idx = find(abs(LME_RFR.Lat - SCal.lat) == min(abs(LME_RFR.Lat - SCal.lat)));
SCal.fCO2 = squeeze(LME_RFR.fCO2(SCal.lon_idx,SCal.lat_idx,:));
SCal.ufCO2 = squeeze(LME_RFR_u.ufCO2(SCal.lon_idx,SCal.lat_idx,:));
SCal.pH = squeeze(LME_RFR.pH(SCal.lon_idx,SCal.lat_idx,:));
SCal.upH = squeeze(LME_RFR_u.upH(SCal.lon_idx,SCal.lat_idx,:));

%% calculate climatology
NGoA.fCO2_clim = nan(12,1);
SGoA.fCO2_clim = nan(12,1);
NCal.fCO2_clim = nan(12,1);
MCal.fCO2_clim = nan(12,1);
SCal.fCO2_clim = nan(12,1);
for m = 1:12
    NGoA.fCO2_clim(m) = mean(NGoA.fCO2(m:12:end));
    SGoA.fCO2_clim(m) = mean(SGoA.fCO2(m:12:end));
    NCal.fCO2_clim(m) = mean(NCal.fCO2(m:12:end));
    MCal.fCO2_clim(m) = mean(MCal.fCO2(m:12:end));
    SCal.fCO2_clim(m) = mean(SCal.fCO2(m:12:end));
end

%% plot fCO2 Time Series
clrs = cbrewer('qual','Set1',5);
f=figure; hold on;
f.Position(3) = 3*f.Position(3);
plot(Time,NGoA.fCO2,'color',clrs(1,:),'linewidth',2);
fill([Time;flipud(Time)],[NGoA.fCO2+NGoA.ufCO2;flipud(NGoA.fCO2-NGoA.ufCO2)],...
    clrs(1,:),'facealpha',0.25,'linestyle','none');
plot(Time,SGoA.fCO2,'color',clrs(2,:),'linewidth',2);
fill([Time;flipud(Time)],[SGoA.fCO2+SGoA.ufCO2;flipud(SGoA.fCO2-SGoA.ufCO2)],...
    clrs(2,:),'facealpha',0.25,'linestyle','none');
plot(Time,NCal.fCO2,'color',clrs(3,:),'linewidth',2);
fill([Time;flipud(Time)],[NCal.fCO2+NCal.ufCO2;flipud(NCal.fCO2-NCal.ufCO2)],...
    clrs(3,:),'facealpha',0.25,'linestyle','none');
plot(Time,MCal.fCO2,'color',clrs(4,:),'linewidth',2);
fill([Time;flipud(Time)],[MCal.fCO2+MCal.ufCO2;flipud(MCal.fCO2-MCal.ufCO2)],...
    clrs(4,:),'facealpha',0.25,'linestyle','none');
plot(Time,SCal.fCO2,'color',clrs(5,:),'linewidth',2);
fill([Time;flipud(Time)],[SCal.fCO2+SCal.ufCO2;flipud(SCal.fCO2-SCal.ufCO2)],...
    clrs(5,:),'facealpha',0.25,'linestyle','none');
datetick('x','yyyy','keeplimits');
exportgraphics(f,'Figures/location_comp.png');
close

%% plot fCO2 Climatology
clrs = cbrewer('qual','Set1',5);
f=figure; hold on;
plot(1:12,NGoA.fCO2_clim-mean(NGoA.fCO2_clim),'color',clrs(1,:),'linewidth',2);
plot(1:12,SGoA.fCO2_clim-mean(SGoA.fCO2_clim),'color',clrs(2,:),'linewidth',2);
plot(1:12,NCal.fCO2_clim-mean(NCal.fCO2_clim),'color',clrs(3,:),'linewidth',2);
plot(1:12,MCal.fCO2_clim-mean(MCal.fCO2_clim),'color',clrs(4,:),'linewidth',2);
plot(1:12,SCal.fCO2_clim-mean(SCal.fCO2_clim),'color',clrs(5,:),'linewidth',2);
plot(0:13,repmat(0,14,1),'k--')
xlim([0 13]);
xticks(1:12);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('{\itf}CO_{2} Anomaly (\muatm)')
exportgraphics(f,'Figures/location_clim_comp.png');
close
max(NGoA.fCO2_clim)-min(NGoA.fCO2_clim)
max(SGoA.fCO2_clim)-min(SGoA.fCO2_clim)
max(NCal.fCO2_clim)-min(NCal.fCO2_clim)
max(MCal.fCO2_clim)-min(MCal.fCO2_clim)
max(SCal.fCO2_clim)-min(SCal.fCO2_clim)

%% plot pH Time Series
clrs = cbrewer('qual','Set1',5);
f=figure; hold on;
f.Position(3) = 3*f.Position(3);
plot(Time,NGoA.pH,'color',clrs(1,:),'linewidth',2);
fill([Time;flipud(Time)],[NGoA.pH+NGoA.upH;flipud(NGoA.pH-NGoA.upH)],...
    clrs(1,:),'facealpha',0.25,'linestyle','none');
plot(Time,SGoA.pH,'color',clrs(2,:),'linewidth',2);
fill([Time;flipud(Time)],[SGoA.pH+SGoA.upH;flipud(SGoA.pH-SGoA.upH)],...
    clrs(2,:),'facealpha',0.25,'linestyle','none');
plot(Time,NCal.pH,'color',clrs(3,:),'linewidth',2);
fill([Time;flipud(Time)],[NCal.pH+NCal.upH;flipud(NCal.pH-NCal.upH)],...
    clrs(3,:),'facealpha',0.25,'linestyle','none');
plot(Time,MCal.pH,'color',clrs(4,:),'linewidth',2);
fill([Time;flipud(Time)],[MCal.pH+MCal.upH;flipud(MCal.pH-MCal.upH)],...
    clrs(4,:),'facealpha',0.25,'linestyle','none');
plot(Time,SCal.pH,'color',clrs(5,:),'linewidth',2);
fill([Time;flipud(Time)],[SCal.pH+SCal.upH;flipud(SCal.pH-SCal.upH)],...
    clrs(5,:),'facealpha',0.25,'linestyle','none');
datetick('x','yyyy','keeplimits');
exportgraphics(f,'Figures/location_comp_pH.png');
close

%% plot averages on map
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([28 62],[194 244]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(parula(18));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itf}CO_{2}';
cbarrow;
% plot land
bordersm('alaska','facecolor',rgb('gray'))
bordersm('continental us','facecolor',rgb('gray'))
bordersm('canada','facecolor',rgb('light grey'))
mlabel off
% plot regions
for n = 1:4
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    z = mean(OAI_grid.(region{n}).fCO2,3,'omitnan')';
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        z,295:10:475,'LineStyle','none');
    alpha 0.5
    clear vars_grid z
end
% plot borders around regions
for n = 1:4
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add locations to plot
scatterm(NGoA.lat,NGoA.lon,100,clrs(1,:),'^','filled');
scatterm(SGoA.lat,SGoA.lon,100,clrs(2,:),'^','filled');
scatterm(NCal.lat,NCal.lon,100,clrs(3,:),'^','filled');
scatterm(MCal.lat,MCal.lon,100,clrs(4,:),'^','filled');
scatterm(SCal.lat,SCal.lon,100,clrs(5,:),'^','filled');
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Location_map.png');

%% plot amplitudes on map
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([28 62],[194 244]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(cmocean('thermal',15));
caxis([1 2.5]);
c.TickLength = 0;
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.String = 'Sea Surface {\itf}CO_{2} Amplitude (\muatm)';
cbarrow;
% plot land
bordersm('alaska','facecolor',rgb('gray'))
bordersm('continental us','facecolor',rgb('gray'))
bordersm('canada','facecolor',rgb('light grey'))
mlabel off
% plot regions
for n = 1:4
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % calculate climatology
    OAI_grid.(region{n}).fCO2_clim = ...
        nan(OAI_grid.(region{n}).dim.x,OAI_grid.(region{n}).dim.y,12);
    for m = 1:12
        OAI_grid.(region{n}).fCO2_clim(:,:,m) = ...
            mean(OAI_grid.(region{n}).fCO2(:,:,m:12:end),3,'omitnan');
    end
    % calculate amplitude
    OAI_grid.(region{n}).fCO2_amp = ...
        max(OAI_grid.(region{n}).fCO2_clim,[],3) - ...
        min(OAI_grid.(region{n}).fCO2_clim,[],3);
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        log10(OAI_grid.(region{n}).fCO2_amp)',1:0.1:2.5,'LineStyle','none');
    alpha 0.5
    clear vars_grid z
end
% plot borders around regions
for n = 1:4
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add locations to plot
scatterm(NGoA.lat,NGoA.lon,100,clrs(1,:),'^','filled');
scatterm(SGoA.lat,SGoA.lon,100,clrs(2,:),'^','filled');
scatterm(NCal.lat,NCal.lon,100,clrs(3,:),'^','filled');
scatterm(MCal.lat,MCal.lon,100,clrs(4,:),'^','filled');
scatterm(SCal.lat,SCal.lon,100,clrs(5,:),'^','filled');
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Location_map.png');

%% plot amplitudes on map (Hawaii)
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([15 25],[190 210]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(cmocean('thermal',15));
caxis([1 2.5]);
c.TickLength = 0;
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.String = 'Sea Surface {\itf}CO_{2} Amplitude (\muatm)';
cbarrow;
% plot land
bordersm('hawaii','facecolor',rgb('gray'))
mlabel off
% plot regions
for n = 11
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    % calculate climatology
    OAI_grid.(region{n}).fCO2_clim = ...
        nan(OAI_grid.(region{n}).dim.x,OAI_grid.(region{n}).dim.y,12);
    for m = 1:12
        OAI_grid.(region{n}).fCO2_clim(:,:,m) = ...
            mean(OAI_grid.(region{n}).fCO2(:,:,m:12:end),3,'omitnan');
    end
    % calculate amplitude
    OAI_grid.(region{n}).fCO2_amp = ...
        max(OAI_grid.(region{n}).fCO2_clim,[],3) - ...
        min(OAI_grid.(region{n}).fCO2_clim,[],3);
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        log10(OAI_grid.(region{n}).fCO2_amp)',1:0.1:2.5,'LineStyle','none');
    alpha 0.5
    clear vars_grid z
end
% plot borders around regions
for n = 11
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% add locations to plot
scatterm(NGoA.lat,NGoA.lon,100,clrs(1,:),'^','filled');
scatterm(SGoA.lat,SGoA.lon,100,clrs(2,:),'^','filled');
scatterm(NCal.lat,NCal.lon,100,clrs(3,:),'^','filled');
scatterm(MCal.lat,MCal.lon,100,clrs(4,:),'^','filled');
scatterm(SCal.lat,SCal.lon,100,clrs(5,:),'^','filled');
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/Location_map.png');

%% print fCO2 means
% NGoA
mod = polyfit(Time,NGoA.fCO2,1);
disp(['NGoA Mean = ' num2str(round(mean(NGoA.fCO2),0))]);
disp(['NGoA Tr. = ' num2str(round(mean(mod(1)*365.25),2))]);
NGoA.fCO2_clim = nan(12,1);
for m = 1:12; NGoA.fCO2_clim(m) = mean(NGoA.fCO2(m:12:end)); end
disp(['NGoA Amp. = ' num2str(round(max(NGoA.fCO2_clim)-min(NGoA.fCO2_clim),0))]);
% SGoA
mod = polyfit(Time,SGoA.fCO2,1);
disp(['SGoA Mean = ' num2str(round(mean(SGoA.fCO2),0))]);
disp(['SGoA Tr. = ' num2str(round(mean(mod(1)*365.25),2))]);
SGoA.fCO2_clim = nan(12,1);
for m = 1:12; SGoA.fCO2_clim(m) = mean(SGoA.fCO2(m:12:end)); end
disp(['SGoA Amp. = ' num2str(round(max(SGoA.fCO2_clim)-min(SGoA.fCO2_clim),0))]);
% NCal
mod = polyfit(Time,NCal.fCO2,1);
disp(['NCal Mean = ' num2str(round(mean(NCal.fCO2),0))]);
disp(['NCal Tr. = ' num2str(round(mean(mod(1)*365.25),2))]);
NCal.fCO2_clim = nan(12,1);
for m = 1:12; NCal.fCO2_clim(m) = mean(NCal.fCO2(m:12:end)); end
disp(['NCal Amp. = ' num2str(round(max(NCal.fCO2_clim)-min(NCal.fCO2_clim),0))]);
% MCal
mod = polyfit(Time,MCal.fCO2,1);
disp(['MCal Mean = ' num2str(round(mean(MCal.fCO2),0))]);
disp(['MCal Tr. = ' num2str(round(mean(mod(1)*365.25),2))]);
MCal.fCO2_clim = nan(12,1);
for m = 1:12; MCal.fCO2_clim(m) = mean(MCal.fCO2(m:12:end)); end
disp(['MCal Amp. = ' num2str(round(max(MCal.fCO2_clim)-min(MCal.fCO2_clim),0))]);
% SCal
mod = polyfit(Time,SCal.fCO2,1);
disp(['SCal Mean = ' num2str(round(mean(SCal.fCO2),0))]);
disp(['SCal Tr. = ' num2str(round(mean(mod(1)*365.25),2))]);
SCal.fCO2_clim = nan(12,1);
for m = 1:12; SCal.fCO2_clim(m) = mean(SCal.fCO2(m:12:end)); end
disp(['SCal Amp. = ' num2str(round(max(SCal.fCO2_clim)-min(SCal.fCO2_clim),0))]);

%% print pH means
% NGoA
mod = polyfit(Time,NGoA.pH,1);
disp(['NGoA Mean = ' num2str(round(mean(NGoA.pH),2))]);
disp(['NGoA Tr. = ' num2str(round(mean(mod(1)*365.25),4))]);
NGoA.pH_clim = nan(12,1);
for m = 1:12; NGoA.pH_clim(m) = mean(NGoA.pH(m:12:end)); end
disp(['NGoA Amp. = ' num2str(round(max(NGoA.pH_clim)-min(NGoA.pH_clim),3))]);
% SGoA
mod = polyfit(Time,SGoA.pH,1);
disp(['SGoA Mean = ' num2str(round(mean(SGoA.pH),2))]);
disp(['SGoA Tr. = ' num2str(round(mean(mod(1)*365.25),4))]);
SGoA.pH_clim = nan(12,1);
for m = 1:12; SGoA.pH_clim(m) = mean(SGoA.pH(m:12:end)); end
disp(['SGoA Amp. = ' num2str(round(max(SGoA.pH_clim)-min(SGoA.pH_clim),3))]);
% NCal
mod = polyfit(Time,NCal.pH,1);
disp(['NCal Mean = ' num2str(round(mean(NCal.pH),2))]);
disp(['NCal Tr. = ' num2str(round(mean(mod(1)*365.25),4))]);
NCal.pH_clim = nan(12,1);
for m = 1:12; NCal.pH_clim(m) = mean(NCal.pH(m:12:end)); end
disp(['NCal Amp. = ' num2str(round(max(NCal.pH_clim)-min(NCal.pH_clim),3))]);
% MCal
mod = polyfit(Time,MCal.pH,1);
disp(['MCal Mean = ' num2str(round(mean(MCal.pH),2))]);
disp(['MCal Tr. = ' num2str(round(mean(mod(1)*365.25),4))]);
MCal.pH_clim = nan(12,1);
for m = 1:12; MCal.pH_clim(m) = mean(MCal.pH(m:12:end)); end
disp(['MCal Amp. = ' num2str(round(max(MCal.pH_clim)-min(MCal.pH_clim),3))]);
% SCal
mod = polyfit(Time,SCal.pH,1);
disp(['SCal Mean = ' num2str(round(mean(SCal.pH),2))]);
disp(['SCal Tr. = ' num2str(round(mean(mod(1)*365.25),4))]);
SCal.pH_clim = nan(12,1);
for m = 1:12; SCal.pH_clim(m) = mean(SCal.pH(m:12:end)); end
disp(['SCal Amp. = ' num2str(round(max(SCal.pH_clim)-min(SCal.pH_clim),3))]);

