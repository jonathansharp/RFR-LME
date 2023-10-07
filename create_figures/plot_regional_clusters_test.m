% Plot groups on map (test)
% define regions
define_regions_eiwg
% set cluster options
set_gmm_options
for n = 1:length(region)
    levs = num_groups(n);
    % load data
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/gridded_clusters_test'],'Clusts_grid');
    figure('visible','on');
    worldmap([min(Preds_grid.(region{n}).lat) max(Preds_grid.(region{n}).lat)],...
        [min(Preds_grid.(region{n}).lon) max(Preds_grid.(region{n}).lon)]);
    pcolorm(Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
        mode(Clusts_grid.(region{n}).groups,3)');
    % pcolorm(Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
    %     mode(Clusts_grid.(region{n}).probabilities.c5,3)');
    colormap(flipud(jet(levs)));
    % colormap([rgb('red');rgb('cyan');rgb('green');rgb('magenta');rgb('blue')]);
    %         colormap(customcolormap([0 1],[rgb('blue'); 1 1 1]));
    % title('Most frequent cluster');
    plot_land('map');
    c=colorbar;
    caxis([0.5 levs+0.5]);
    c.Label.String = 'Cluster';
    % c.Label.String = 'C5 Probability';
    c.TickLength = 0;
    c.Ticks = 1:levs;
    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/' region{n} '_Clusts_gmm_var_' num2str(levs) '.png']);
    close all
end
% clean up
clear
