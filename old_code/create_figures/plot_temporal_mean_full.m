% Plot temporal mean of spatial surface variable over entire region
% 
% Written by J.D. Sharp: 10/18/22
% Last updated by J.D. Sharp: 1/18/23
% 

function plot_temporal_mean_full(zmin,zmax,clrmp,varname,lab,region,lme_shape,lme_idx)

% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south','ffacecolor','w');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','south','Position',[0.45 0.2 0.3 0.025]);
caxis([zmin zmax]);
colormap(clrmp);
c.FontSize = 10;
c.FontWeight = 'bold';
c.TickLength = 0;
c.Label.String = lab;
cbarrow;
% plot regions
for n = 1:length(region)
    if any(strcmp(varname,{'DIC' 'fCO2' 'pCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' ...
            'CO3' 'RF' 'ufCO2' 'upCO2' 'upH' 'uOmA' 'uOmC' 'uRF' 'uTA' 'TA_DIC'}))
        type = 'OAI_grid';
        vars_grid = load(['Data/' region{n} '/ML_fCO2'],type);
        if strcmp(varname,'H')
            vars_grid.(type).(region{n}).(varname) = ...
                (10^9).*vars_grid.(type).(region{n}).(varname);
        end
    else
        type = 'Preds_grid';
        vars_grid = load(['Data/' region{n} '/gridded_predictors'],type);
    end
    if ~strcmp(varname,'TA_DIC')
        z = mean(vars_grid.(type).(region{n}).(varname),3,'omitnan')';
    else
        z = mean(vars_grid.(type).(region{n}).TA./ ...
            vars_grid.(type).(region{n}).DIC,3,'omitnan')';
    end
    contourfm(vars_grid.(type).(region{n}).lat,vars_grid.(type).(region{n}).lon,...
        z,zmin:(zmax-zmin)/200:zmax,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
plot_lme_borders(region,lme_shape,lme_idx);
% plot land
plot_land('map');
%mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
% exportgraphics(gcf,['Figures/full/' varname '.png']);
export_fig(gcf,['Figures/full/' varname '.png'],'-transparent');
close
