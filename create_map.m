% function to create animations
function create_map(var,type,lat,lon,z,cmap,colorlims,colorlabel,units,paths)

    % establish figure
    f = figure; set(f,'color','w','visible','off');

    % plot map
    m_proj('Robinson','lat',[min(lat) max(lat)],'lon',[min(lon) max(lon)]);
    m_pcolor(lon,lat,z');
    colormap(cmap);
    m_coast('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    clim(colorlims);
    c=colorbar;
    c.Limits = colorlims;
    c.Label.String = [colorlabel ' (' units ')'];
    c.TickLength = 0;
    % save frame
    if ~isfolder([paths.figures '/' var '/' type]); mkdir([paths.figures '/' var '/' type]); end
    export_fig(f,[paths.figures '/' var '/' type '/' strrep(colorlabel,' ','_') '.png'],'-transparent');


end