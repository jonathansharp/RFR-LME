% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

for n = 1:length(region)
    % load data
    load(['Data/' region{n} '/gridded_pco2_no_moorings'],'SOCAT_grid');
    % plot the percentage of grid cells with data in each region
    figure('visible','on');
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
end
