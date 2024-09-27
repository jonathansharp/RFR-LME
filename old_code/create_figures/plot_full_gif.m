% Plot temporal mean of gridded pCO2 over entire region
% 
% Written by J.D. Sharp: 10/28/22
% Last updated by J.D. Sharp: 1/18/23
% 

function plot_full_gif(zmin,zmax,clrmp,varname,lab,region,lme_shape,lme_idx)

% initialize figure
figure('visible','off');
box on; hold on;
set(gcf,'position',[100 100 900 600],'color','white');
% initialize axis
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gca,'fontsize',16);
% colorbar
c=colorbar('location','southoutside');
c.TickLength = 0;
c.Label.String = lab;
colormap(clrmp);
caxis([zmin,zmax]);
cbarrow;
% define filename
filename = ['Figures/full/' varname '_monthly.gif'];
f=1;
% loop through months
for m = 1:300
    % initialize axis
    worldmap([-18 82],[140 302]);
    setm(gca,'MapProjection','robinson','MLabelParallel','south');
    set(gca,'fontsize',16);    
    % plot regions
    for n = 1:length(region)
        if any(strcmp(varname,{'DIC' 'fCO2' 'ufCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF' 'pCO2'}))
            type = 'OAI_grid';
            vars_grid = load(['Data/' region{n} '/ML_fCO2'],type);
        elseif any(strcmp(varname,{'fco2_ave_wtd'}))
            type = 'SOCAT_grid';
            vars_grid = load(['Data/' region{n} '/gridded_pco2'],type);
        else
            type = 'Preds_grid';
            vars_grid = load(['Data/' region{n} '/gridded_predictors'],type);
        end
        z = vars_grid.(type).(region{n}).(varname)(:,:,m)';
%         contourfm(vars_grid.(type).(region{n}).lat,...
%                   vars_grid.(type).(region{n}).lon,...
%                   z,zmin:(zmax-zmin)/200:zmax,'LineStyle','none');
        pcolorm(vars_grid.(type).(region{n}).lat,...
                  vars_grid.(type).(region{n}).lon,z);
    end
    % plot borders around regions
    plot_lme_borders(region,lme_shape,lme_idx);
    % plot land
    plot_land('map');
    mlabel off
    % add text
    textm(38,152,[num2str(vars_grid.(type).(region{n}).month_of_year(m)) '/' ...
        num2str(vars_grid.(type).(region{n}).year(m))],'fontsize',16);
    % clean up
    clear vars_grid z
    % save frame
    drawnow;
    frame = getframe(gcf);
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