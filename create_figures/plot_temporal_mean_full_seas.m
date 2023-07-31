% Plot temporal mean of spatial surface variable over entire region (seasonal)
% 
% Written by J.D. Sharp: 1/23/23
% Last updated by J.D. Sharp: 1/23/23
% 

function plot_temporal_mean_full_seas(zmin,zmax,clrmp,varname,lab,region,lme_shape,lme_idx)

seas_lab = {'DJF' 'MAM' 'JJA' 'SON'};
seas = [1,2,12;3,4,5;6,7,8;9,10,11];

for ms = 1:4
    
    % initialize figure
    figure('visible','off'); box on; hold on;
    worldmap([-18 82],[140 302]);
    setm(gca,'MapProjection','robinson','MLabelParallel','south');
    set(gcf,'position',[100 100 900 600]);
    set(gca,'fontsize',16);
    % figure properties
    c=colorbar('location','southoutside');
    colormap(clrmp);
    caxis([zmin,zmax]);
    c.TickLength = 0;
    c.Label.String = lab;
    cbarrow;
    % plot regions
    for n = 1:length(region)
        if any(strcmp(varname,{'DIC' 'fCO2' 'pCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'}))
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
        z = mean(cat(3,vars_grid.(type).(region{n}).(varname)(:,:,seas(ms,1):12:end),...
                       vars_grid.(type).(region{n}).(varname)(:,:,seas(ms,2):12:end),...
                       vars_grid.(type).(region{n}).(varname)(:,:,seas(ms,3):12:end)),3,'omitnan')';
        contourfm(vars_grid.(type).(region{n}).lat,vars_grid.(type).(region{n}).lon,...
            z,zmin:(zmax-zmin)/200:zmax,'LineStyle','none');
        clear vars_grid z
    end
    % plot borders around regions
    for n = 1:length(region)
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
    mlabel off
    % save figure
    if ~isfolder('Figures/full'); mkdir('Figures/full'); end
    exportgraphics(gcf,['Figures/full/' varname '_' seas_lab{ms} '.png']);
    close

end