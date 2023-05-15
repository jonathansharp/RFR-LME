% Plot temporal mean of spatial surface variable
% 
% Written by J.D. Sharp: 8/24/22
% Last updated by J.D. Sharp: 9/15/22
% 

function plot_temporal_mean(lim,dim,lat,lon,z,clrmp,varname,lab,reg)

zmean = mean(z,3,'omitnan');
figure('visible','off'); worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
if ~contains(varname,'_var') & ~startsWith(varname,'u')
    zmax = mean(zmean(:),'omitnan')+3.*std(zmean(:),[],'omitnan');
    zmin = mean(zmean(:),'omitnan')-3.*std(zmean(:),[],'omitnan');
else
    zmax = mean(zmean(:),'omitnan')+5.*std(zmean(:),[],'omitnan');
    zmin = 0;
end
if ~isnan(zmin) && ~isnan(zmax)
    pcolorm(repmat(lat',dim.x,1),repmat(lon,1,dim.y),...
            mean(z,3,'omitnan'),20);
end
plot_land('map');
c=colorbar;
colormap(clrmp);
if ~isnan(zmin) && ~isnan(zmax)
    if zmin == 0 && zmax == 0
        caxis([0 1]);
    else
        caxis([zmin,zmax]);
    end
end
c.TickLength = 0;
c.Label.String = lab;
cbarrow;

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,['Figures/' reg '_' varname '.png']);
close
