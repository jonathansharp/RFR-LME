% Plot temporal mean of
% 
% Written by J.D. Sharp: 11/30/22
% Last updated by J.D. Sharp: 1/13/22
% 

function plot_temporal_mean_delfco2_full(zmin,zmax,lev,varname,lab,region,lme_shape,lme_idx)

figure('visible','off'); box on; hold on;
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
xlim([169 305]); ylim([5 82]);
xlabel('');
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    z = mean(Preds_grid.(region{n}).(varname),3,'omitnan')';
    h=imagesc(Preds_grid.(region{n}).lon,Preds_grid.(region{n}).lat,z);
    set(h,'alphadata',~isnan(z));
    clear Preds_grid
end
for n = 1:length(region)
    for r = 1:length(lme_idx.(region{n}))
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})(r)).X);
        plot(tmp_lon,lme_shape(lme_idx.(region{n})(r)).Y,'k','linewidth',1);
    end
end
plot_land('xy');
c=colorbar;
colormap(parula(lev));
caxis([zmin,zmax]);
c.TickLength = 0;
c.Label.String = lab;
cbarrow;

% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,['Figures/full/' varname '.png']);
close