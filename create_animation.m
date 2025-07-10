% function to create animations
function create_animation(var,type,time,lat,lon,z,cmap,colorlims,colorlabel,units)

    % establish figure
    f = figure; set(f,'color','w','visible','off');

    % loop through months and plot each one
    for m = 1:length(time)
        m_proj('Robinson','lat',[min(lat) max(lat)],'lon',[min(lon) max(lon)]);
        m_pcolor(lon,lat,z(:,:,m)');
        title(gca,extractAfter(datestr(time(m)),'-'));
        colormap(cmap);
        m_coast('patch',[0.7 0.7 0.7]);
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[]);
        clim(colorlims);
        c=colorbar;
        c.Limits = colorlims;
        c.Label.String = [colorlabel ' (' units ')'];
        c.TickLength = 0;
        % save frame
        if ~isfolder(['Figures/' var '/' type '/Monthly']); mkdir(['Figures/' var '/' type '/Monthly']); end
        exportgraphics(f,['Figures/' var '/' type '/Monthly/m' num2str(m) '_' strrep(colorlabel,' ','_') '.png']);
        % capture frame
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if m == 1
            imwrite(imind,cm,['Figures/' var '/' type '/' ...
                strrep(colorlabel,' ','_') '_Animation.gif'],...
                'gif','Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind,cm,['Figures/' var '/' type '/' ...
                strrep(colorlabel,' ','_') '_Animation.gif'],...
                'gif','WriteMode','append','DelayTime',0.1);
        end
        clf; % clear frame
    end

end
