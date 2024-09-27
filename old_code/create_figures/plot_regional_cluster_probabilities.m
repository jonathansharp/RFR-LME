% Plot probability for each group on map
% define regions
define_regions_eiwg
% set cluster options
set_gmm_options
for n = 1:length(region)
    for g = 1:num_groups(n)
        levs = num_groups(n);
        % load data
        load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
        load(['Data/' region{n} '/gridded_clusters'],'Clusts_grid');
        figure('visible','on');
        worldmap([min(Preds_grid.(region{n}).lat) max(Preds_grid.(region{n}).lat)],...
            [min(Preds_grid.(region{n}).lon) max(Preds_grid.(region{n}).lon)]);
        pcolorm(Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
            mode(Clusts_grid.(region{n}).probabilities.(['c' num2str(g)]),3)');
        clrs = jet(levs);
        colormap(customcolormap([0 1],[clrs(g,:); 1 1 1]));
        plot_land('map');
        c=colorbar;
        c.Label.String = ['C' num2str(g) ' Probability'];
        c.TickLength = 0;
        c.Ticks = 1:levs;
        % save figure
        if ~isfolder('Figures'); mkdir('Figures'); end
        exportgraphics(gcf,['Figures/' region{n} '_Clusts_gmm_var_' ...
            num2str(g) '_probs.png']);
        close all
    end
end
% clean up
clear
