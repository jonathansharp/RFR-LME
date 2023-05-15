% Load Variables
% 
% This script loads satellite and reanalysis variables for 18 US large
% Marine Ecosystems and regrids them in preparation for clustering and
% machine learning algorithm fits.
% 
% Written by J.D. Sharp: 8/19/22
% Last updated by J.D. Sharp: 5/2/23
% 

%% define path
filepath = pwd;

%% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');

%% display status
disp('Loading predictor variables');

%% Duplicate SOCAT grid for predictors grid
Preds_grid.lon = SOCAT_grid.lon;
Preds_grid.lat = SOCAT_grid.lat;
Preds_grid.lim = SOCAT_grid.lim;
Preds_grid.dim = SOCAT_grid.dim;
Preds_grid.month = SOCAT_grid.month;
Preds_grid.year = SOCAT_grid.year;
Preds_grid.month_of_year = SOCAT_grid.month_of_year;
clear SOCAT_grid

%% Calculate distance from coast
Preds_grid.Dist = ...
    dist2coast(repmat(Preds_grid.lat',Preds_grid.dim.x,1),...
    repmat(Preds_grid.lon,1,Preds_grid.dim.y));

%% Obtain sea surface salinity from GLORYS
Preds_grid.SSS = ...
    import_SSS_GLORYS_all(Preds_grid.lat,Preds_grid.lon,...
                    Preds_grid.month,filepath);
% test plot
figure;
worldmap([Preds_grid.lim.latmin Preds_grid.lim.latmax],...
    [Preds_grid.lim.lonmin Preds_grid.lim.lonmax]);
pcolorm(Preds_grid.lat,Preds_grid.lon,mean(Preds_grid.SSS,3,'omitnan')');
plot_land('map');
colorbar;

%% Obtain sea surface height from CMEMS satellite product
Preds_grid.SSH = ...
    import_SSH_CMEMS_all(Preds_grid.lat,Preds_grid.lon,filepath);
% test plot
figure;
worldmap([Preds_grid.lim.latmin Preds_grid.lim.latmax],...
    [Preds_grid.lim.lonmin Preds_grid.lim.lonmax]);
pcolorm(Preds_grid.lat,Preds_grid.lon,mean(Preds_grid.SSH,3,'omitnan')');
plot_land('map');
colorbar;

%% Obtain sea surface temperature from OISSTv2
Preds_grid.SST = ...
    import_SST_OISST_all(Preds_grid.lat,Preds_grid.lon,filepath);
% test plot

%% Obtain sea surface ice concentration from OISSTv2
Preds_grid.IceC = ...
    import_Ice_OISST(Preds_grid.lat,Preds_grid.lon,filepath);
% test plot

%% Obtain sea surface chlorophyll from satellite measurements
Preds_grid.CHL = ...
    import_CHL_NASA(Preds_grid.lat,...
                    Preds_grid.lon,...
                    Preds_grid.month,...
                    Preds_grid.idxspc(:,:,1),filepath);
Preds_grid.CHL(Preds_grid.CHL<0) = 0.0001;
% test plot

%% Obtain wind speed from ERA5 re-analysis
Preds_grid.WindSpeed = ...
    import_Winds_ERA5(Preds_grid.lat,Preds_grid.lon,Preds_grid.month,filepath);
% test plot

%% Obtain bathymetry from ETOPO2
Preds_grid.Bathy = ...
    import_Bathy_ETOPO2(Preds_grid.lat,Preds_grid.lon,filepath);
% negative trap for bathymetry
Preds_grid.Bathy(Preds_grid.Bathy < 0) = 0;
% test plot

%% Obtain mixed layer depth
Preds_grid.MLD = ...
    import_MLD_GLORYS(Preds_grid.lat,Preds_grid.lon,...
                     Preds_grid.month,filepath);
% Negative trap for MLD
Preds_grid.MLD(Preds_grid.MLD<0) = 10;
% test plot

%% Obtain atmospheric pressure from NCEP
Preds_grid.mslp = ...
    import_mslp_NCEP(Preds_grid.lat,Preds_grid.lon,...
                     Preds_grid.month,filepath);
% test plot

%% Obtain atmospheric pCO2 from NOAA MBL product
Preds_grid.pCO2_atm = ...
    import_pCO2_MBL(Preds_grid.lat,Preds_grid.lon,...
                    Preds_grid.month,Preds_grid.SSS,...
                    Preds_grid.SST,Preds_grid.mslp,filepath);
% test plot

%% Save gridded predictor data
save('Data/gridded_predictors','Preds_grid','-v7.3');

% clean up
clear Preds_grid SOCAT_grid

end

% clean up
clear filepath

%% Plot number of all grid cells with data
% load SOCAT grid
load('Data/socat_gridded','SOCAT_grid');
% initialize figure
figure('visible','off'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
mycolormap = jet(21);
mycolormap(1,:) = 1;
colormap(mycolormap);
caxis([-1 41]);
c.TickLength = 0;
c.Label.String = 'Total Months Represented';
cbarrow('right');
% plot background
pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,SOCAT_grid.num_months');
alpha 0.3
% clear SOCAT grid
clear SOCAT_grid
% plot regions
for n = 1:length(region)
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    pcolorm(SOCAT_grid.lat,SOCAT_grid.lon,...
        SOCAT_grid.num_months');
    clear SOCAT_grid
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
exportgraphics(gcf,'Figures/full/num_obs.png');
close
% clean up
clear n h r tmp_lon c mycolormap 

%% plot surface variables across full region
plot_temporal_mean_full(29.75,38.25,0.5,cmocean('haline',17),'SSS','Sea Surface Salinity',region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.015,0.105,0.01,parula(12),'SSH','Sea Surface Height Anomaly (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(1,31,2,cmocean('thermal',15),'SST',['Sea Surface Temperature (' char(176) 'C)'],region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.025,1.025,0.05,cmocean('tempo',21),'IceC','Sea Ice Concentration Fraction',region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.025,1.025,0.05,cmocean('algae',21),'CHL','Sea Surface Chlorophyll (log_{10})',region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.5,12.5,1,cmocean('amp',13),'WindSpeed','Wind Speed (m s^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(-125,6125,250,cmocean('deep',25),'Bathy','Bottom Depth (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(-2.5,52.5,5,jet(11),'MLD','Mixed Layer Depth (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(0.990,1.010,0.001,cmocean('dense',21),'mslp','Sea Level Pressure (atm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(369.5,390.5,1,cmocean('solar',21),'pCO2_atm','Atmospheric pCO_{2} (\muatm)',region,lme_shape,lme_idx)

%% plot gifs of surface variables across full region
plot_full_gif(29.75,38.25,cmocean('haline',17),'SSS','Sea Surface Salinity',region,lme_shape,lme_idx)
plot_full_gif(-0.015,0.105,parula(12),'SSH','Sea Surface Height Anomaly (m)',region,lme_shape,lme_idx)
plot_full_gif(1,31,cmocean('thermal',15),'SST',['Sea Surface Temperature (' char(176) 'C)'],region,lme_shape,lme_idx)
plot_full_gif(-0.025,1.025,cmocean('tempo',21),'IceC','Sea Ice Concentration Fraction',region,lme_shape,lme_idx)
plot_full_gif(-0.025,1.025,cmocean('algae',21),'CHL','Sea Surface Chlorophyll (log_{10})',region,lme_shape,lme_idx)
plot_full_gif(-0.5,12.5,cmocean('amp',13),'WindSpeed','Wind Speed (m s^{-1})',region,lme_shape,lme_idx)
plot_full_gif(-2.5,52.5,jet(11),'MLD','Mixed Layer Depth (m)',region,lme_shape,lme_idx)
plot_full_gif(0.990,1.010,cmocean('dense',21),'mslp','Sea Level Pressure (atm)',region,lme_shape,lme_idx)
plot_full_gif(369.5,390.5,cmocean('solar',21),'pCO2_atm','Atmospheric pCO_{2} (\muatm)',region,lme_shape,lme_idx)
