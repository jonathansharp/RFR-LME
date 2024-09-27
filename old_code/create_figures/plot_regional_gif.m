% Plot temporal mean of gridded pCO2 over entire region
% 
% Written by J.D. Sharp: 10/28/22
% Last updated by J.D. Sharp: 1/18/23
% 

function plot_regional_gif(lim,lat,lon,z,clrmp,varname,lab,year,...
    month_of_year,reg,shape)

% initialize figure
h=figure('visible','off');
set(h,'color','white');
box on; hold on;
% determine mean and limits
zmean = mean(z,3,'omitnan');
if ~contains(varname,'_var') && ~startsWith(varname,'u')
    zmax = mean(zmean(:),'omitnan')+2.*std(zmean(:),[],'omitnan');
    zmin = mean(zmean(:),'omitnan')-2.*std(zmean(:),[],'omitnan');
else
    zmax = mean(zmean(:),'omitnan')+4.*std(zmean(:),[],'omitnan');
    zmin = 0;
end
% colorbar properties
c=colorbar('location','eastoutside');
c.TickLength = 0;
c.Label.String = lab;
colormap(clrmp);
caxis([zmin,zmax]);
cbarrow;
% define filename
filename = ['Figures/' reg '_' varname '_monthly.gif'];
f=1;
% loop through months
for m = 1:length(month_of_year)
    % initialize axis
    worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
    % plot
    %pcolorm(lat,lon,z(:,:,m)');
    contourfm(lat,lon,z(:,:,m)',zmin:(zmax-zmin)/200:zmax,'LineStyle','none');
    % plot land
    plot_land('map');
    % plot border
    plotm(shape.Y,shape.X,'k','linewidth',1);
    % add text
    title([num2str(month_of_year(m)) '/' num2str(year(m))],'fontsize',16);
    % save frame
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if f == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    cla;
    f=f+1;
end

close