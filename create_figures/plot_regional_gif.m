% Plot temporal mean of gridded pCO2 over entire region
% 
% Written by J.D. Sharp: 10/28/22
% Last updated by J.D. Sharp: 1/18/23
% 

function plot_regional_gif(lim,lat,lon,z,clrmp,varname,lab,year,month_of_year,reg)

% initialize figure
h=figure('visible','on');
set(h,'color','white');
box on; hold on;
worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
% determine mean and limits
zmean = mean(z,3,'omitnan');
if ~contains(varname,'_var')
    zmax = mean(zmean(:),'omitnan')+3.*std(zmean(:),[],'omitnan');
    zmin = mean(zmean(:),'omitnan')-3.*std(zmean(:),[],'omitnan');
else
    zmax = mean(zmean(:),'omitnan')+5.*std(zmean(:),[],'omitnan');
    zmin = 0;
end
% colorbar properties
c=colorbar('location','eastoutside');
c.TickLength = 0;
c.Label.String = lab;
cbarrow;
colormap(clrmp);
caxis([zmin,zmax]);
% define filename
filename = ['Figures/' reg '_' varname '_monthly.gif'];
f=1;
% loop through months
for m = 1:288
    % initialize axis
%     worldmap([lim.latmin lim.latmax],[lim.lonmin lim.lonmax]);
%     setm(gca,'MapProjection','robinson');
%     set(gca,'fontsize',16);
    % plot
    p=pcolorm(lat,lon,z(:,:,m)');
    % plot land
    plot_land('map');
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
    f=f+1;
end

close