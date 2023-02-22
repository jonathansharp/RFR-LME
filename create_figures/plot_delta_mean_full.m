% Plot temporal mean of spatial surface variable over entire region
% 
% Written by J.D. Sharp: 10/18/22
% Last updated by J.D. Sharp: 1/31/23
% 

function plot_delta_mean_full(zmin,zmax,cmap_type,cmap_name,cmap_segs,zero_piv,num_groups,varname,lab,region,lme_shape,lme_idx)

% initialize figure
figure('visible','off'); hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar;
caxis([zmin zmax]);
if strcmp(cmap_type,'cmocean')
    if zero_piv == 1
        colormap(cmocean(cmap_name,cmap_segs,'pivot',0));
    else
        colormap(cmocean(cmap_name,cmap_segs));
    end
elseif strcmp(cmap_type,'stnd')
    colormap(cmap_name);
end
c.TickLength = 0;
c.Label.String = lab;
cbarrow;
% plot regions
type = 'Val';
for n = 1:length(region)
for en = 1:size(num_groups,2)    
    vars_grid_temp = load(['Data/' region{n} '/us_lme_model_evals_' num2str(en)],type);
    vars_grid.([type num2str(en)]).(region{n}) = ...
        vars_grid_temp.(type).(region{n});
    clear vars_grid_temp
end
vars_grid.(type).(region{n}).(varname) = ...
    mean(cat(4,vars_grid.Val1.(region{n}).(varname)),4);
z = mean(vars_grid.(type).(region{n}).(varname),3,'omitnan')';
lat = vars_grid.([type num2str(en)]).(region{n}).lat;
lon = vars_grid.([type num2str(en)]).(region{n}).lon;
pcolorm(lat,lon,z)
% plot borders around regions
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% plot land
plot_land('map');
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,['Figures/full/' varname '.png']);
close