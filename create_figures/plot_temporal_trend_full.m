% Plot temporal trend of spatial surface variable over entire region
% 
% Written by J.D. Sharp: 10/18/22
% Last updated by J.D. Sharp: 1/18/23
% 

function plot_temporal_trend_full(zmin,zmax,clrmp,varname,lab,region,lme_shape,lme_idx)

% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','south','Position',[0.45 0.2 0.3 0.025]);
caxis([zmin zmax]);
colormap(cmocean('balance','pivot',0));
c.FontSize = 10;
c.FontWeight = 'bold';
c.TickLength = 0;
if varname == 'pH'
    c.TickLabels = {'-2' '-1' '0' '1' '2'};
end
c.Label.String = lab;
cbarrow;
% plot regions
for n = 1:length(region)
    if any(strcmp(varname,{'DIC' 'fCO2' 'pCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF' 'ufCO2' 'upCO2' 'TA_DIC'}))
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
    tr = nan(vars_grid.(type).(region{n}).dim.x,vars_grid.(type).(region{n}).dim.y);
    for a = 1:vars_grid.(type).(region{n}).dim.x
        for b = 1:vars_grid.(type).(region{n}).dim.y
            if sum(~isnan(squeeze(vars_grid.(type).(region{n}).(varname)(a,b,:))))>1 % if more than 200 months with observations
            [~,~,x,~] = leastsq2(vars_grid.(type).(region{n}).month,...
                squeeze(vars_grid.(type).(region{n}).(varname)(a,b,:)),0,2,[6 12]);
            tr(a,b) = x(2)*12;
            else
                tr(a,b) = NaN;
            end
        end
    end
    contourfm(vars_grid.(type).(region{n}).lat,vars_grid.(type).(region{n}).lon,...
        tr',zmin:(zmax-zmin)/200:zmax,'LineStyle','none');
    clear vars_grid tr
end
% plot borders around regions
plot_lme_borders(region,lme_shape,lme_idx);
% plot land
plot_land('map');
mlabel off
% save figure
if ~isfolder('Figures/full'); mkdir('Figures/full'); end
exportgraphics(gcf,['Figures/full/' varname '_tr.png']);
close
