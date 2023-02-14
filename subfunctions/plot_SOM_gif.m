n=4;
load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');

%% Plot cluster GIF
h=figure('visible','off');
set(h,'color','white');
axis tight manual
filename = ['Figures/' region{n} '/SOM_monthly.gif'];
f=1;
for m = 1:Preds_grid.(region{n}).dim.z
    % plot clusters
    worldmap([min(Preds_grid.(region{n}).lat) max(Preds_grid.(region{n}).lat)],...
        [min(Preds_grid.(region{n}).lon) max(Preds_grid.(region{n}).lon)]);
    pcolorm(Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
        Preds_grid.(region{n}).groups(:,:,m)');
    colormap(jet(9));
    title([num2str(Preds_grid.(region{n}).month_of_year(m)) '/' ...
        num2str(Preds_grid.(region{n}).year(m))]);
    plot_land('map');
    c=colorbar;
    caxis([0.5 9.5]);
    c.Label.String = 'Cluster';
    c.TickLength = 0;
    % save frame
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if f == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    f=f+1;
end

% save figure
% if ~isfolder(['Figures/' region{n}]); mkdir(['Figures/' region{n}]); end
% exportgraphics(gcf,['Figures/' region{n} '/SOM.png']);
% close all