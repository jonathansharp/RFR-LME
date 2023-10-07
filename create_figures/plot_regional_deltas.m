% Plot gridded delta pCO2 for RFRs

% define regions
define_regions_eiwg

for n = 1:length(region)
    % load data
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/us_lme_model_evals'],'Val');
    % create plot
    figure('visible','off');
    worldmap([Preds_grid.(region{n}).lim.latmin ...
        Preds_grid.(region{n}).lim.latmax],...
       [Preds_grid.(region{n}).lim.lonmin ...
        Preds_grid.(region{n}).lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(Preds_grid.(region{n}).lat',Preds_grid.(region{n}).dim.x,1),...
            repmat(Preds_grid.(region{n}).lon,1,Preds_grid.(region{n}).dim.y),...
            mean(abs(Val.(region{n}).delta_rfr_grid),3,'omitnan'));
    plot_land('map');
    c=colorbar;
    colormap(cmocean('tempo',16));
    caxis([-1 31]);
    c.TickLength = 0;
    c.Label.String = '\Delta{\itf}CO_{2} (\muatm)';
    cbarrow('up');
    clear c
    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/' region{n} '_delta_fCO2_RFR_gridded.png']);
    close
end
% clean up
clear