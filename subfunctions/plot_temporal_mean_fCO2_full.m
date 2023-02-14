% Plot temporal mean of gridded fCO2 over entire region
% 
% Written by J.D. Sharp: 10/28/22
% Last updated by J.D. Sharp: 1/16/23
% 

function plot_temporal_mean_fCO2_full(zmin,zmax,zsp,lev,varname,lab,region,lme_shape,lme_idx)

% initialize figure
figure('visible','off');
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
box on; hold on;
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar;
colormap(parula(lev));
caxis([zmin,zmax]);
c.TickLength = 0;
c.Label.String = lab;
cbarrow;
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    z = mean(OAI_grid.(region{n}).(varname),3,'omitnan')';
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        z,zmin:zsp:zmax,'LineStyle','none');
    clear OAI_grid
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
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,['Figures/full/' varname '.png']);
close