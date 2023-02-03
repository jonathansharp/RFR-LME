% Predict TA on grid and perform CO2 system calculations
% 
% This script uses ESPER algorithms to estimate surface TA on a grid for US
% large Marine Ecosystems.
% 
% Written by J.D. Sharp: 1/19/23
% Last updated by J.D. Sharp: 2/1/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions

for n = 1:length(region)

    %% display status
    disp(['Calculating CO2 System (' region{n} ')']);

    %% load gridded fCO2, predictors, and models
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');

    %% predict TA (and nutrients) using ESPER
    % process ESPER predictors
    lon = repmat(Preds_grid.(region{n}).lon,1,...
        Preds_grid.(region{n}).dim.y,Preds_grid.(region{n}).dim.z);
    lat = repmat(Preds_grid.(region{n}).lat',...
        Preds_grid.(region{n}).dim.x,1,Preds_grid.(region{n}).dim.z);
    lon = lon(Preds_grid.(region{n}).idxspc);
    lat = lat(Preds_grid.(region{n}).idxspc);
    depth = ones(size(lon));
    sal = Preds_grid.(region{n}).SSS(Preds_grid.(region{n}).idxspc);
    tmp = Preds_grid.(region{n}).SST(Preds_grid.(region{n}).idxspc);
    % predict TA using ESPER-Mixed
    [data,u_data] = ESPER_Mixed([1 4 6],[lon lat depth],[sal tmp],[1 2],'Equations',8);
    % pre-allocate TA and uTA
    OAI_grid.(region{n}).TA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uTA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).P = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).Si = nan(size(Preds_grid.(region{n}).idxspc));
    % add TA and uTA to OA grid
    OAI_grid.(region{n}).TA(Preds_grid.(region{n}).idxspc) = data.TA;
    OAI_grid.(region{n}).uTA(Preds_grid.(region{n}).idxspc) = u_data.TA;
    OAI_grid.(region{n}).P(Preds_grid.(region{n}).idxspc) = data.phosphate;
    OAI_grid.(region{n}).Si(Preds_grid.(region{n}).idxspc) = data.silicate;

    %% plot estimated TA
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).TA,...
        cmocean('haline',24),'TA',...
        'Surface {\itA}_{T} (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).P,...
        cmocean('haline',24),'P',...
        'Surface Phosphate (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).Si,...
        cmocean('haline',24),'Si',...
        'Surface Silicate (\mumol kg^{-1})',region{n});

    %% calculate other OA Indicators
    %fCO2 = sal = OAI_grid.(region{n}).fCO2(OAI_grid.(region{n}).idxspc);
    fCO2 = OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).idxspc);
    carb_system = CO2SYS(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,1,10,1,2,2);
    % pre-allocate OA indicators
    OAI_grid.(region{n}).DIC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).pH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).OmA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).OmC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).H = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).CO3 = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).RF = nan(size(Preds_grid.(region{n}).idxspc));
    % add OA indicators to OA grid
    OAI_grid.(region{n}).DIC(Preds_grid.(region{n}).idxspc) = carb_system(:,2);
    OAI_grid.(region{n}).pH(Preds_grid.(region{n}).idxspc) = carb_system(:,3);
    OAI_grid.(region{n}).OmA(Preds_grid.(region{n}).idxspc) = carb_system(:,18);
    OAI_grid.(region{n}).OmC(Preds_grid.(region{n}).idxspc) = carb_system(:,17);
    OAI_grid.(region{n}).H(Preds_grid.(region{n}).idxspc) = 10.^-carb_system(:,3);
    OAI_grid.(region{n}).CO3(Preds_grid.(region{n}).idxspc) = carb_system(:,7);
    OAI_grid.(region{n}).RF(Preds_grid.(region{n}).idxspc) = carb_system(:,16);
    
    %% plot estimated OA indicators
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).DIC,...
        flipud(jet(20)),'DIC','Dissolved Inorganic Carbon',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).pH,...
        flipud(jet(20)),'pH','Surface pH_{T}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).OmA,...
        flipud(jet(16)),'OmA','Surface \Omega_{A}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).OmC,...
        flipud(jet(16)),'OmC','Surface \Omega_{C}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,(10^9).*OAI_grid.(region{n}).H,...
        jet(20),'H','Surface [H^{+}] (nmol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).CO3,...
        jet(15),'CO3','Surface [CO_{3}^{2-}] (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).RF,...
        flipud(jet(28)),'RF','Surface RF',region{n});

    %% save estimated OA grid
    save(['Data/' region{n} '/ML_fCO2'],'OAI_grid','-v7.3');

    %% clean up
    clear carb_system depth fCO2 lat lon OAI_grid Preds_grid sal TA tmp uTA

end

%% plot time series
OA_time_series
OA_summary_stats

%% plot OA indicators across full region
plot_temporal_mean_full(1900,2500,25,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(1700,2300,25,flipud(jet(15)),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
plot_temporal_mean_full(7.9,8.2,0.02,jet(15),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,5,0.25,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
plot_temporal_mean_full(1,6,0.25,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
plot_temporal_mean_full(5,10,0.25,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
plot_temporal_mean_full(50,250,10,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
plot_temporal_mean_full(8,15,0.5,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)

%% plot gifs of OA indicators across full region
plot_full_gif(1900,2500,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_full_gif(1700,2300,jet(15),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
plot_full_gif(7.9,8.2,flipud(jet(15)),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
plot_full_gif(0,5,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
plot_full_gif(1,6,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
plot_full_gif(5,10,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
plot_full_gif(50,250,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
plot_full_gif(8,16,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)

%% plot OA indicators across full region (seasonally)
plot_temporal_mean_full_seas(1900,2500,25,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(1700,2300,25,flipud(jet(15)),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(7.9,8.2,0.02,flipud(jet(15)),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(0,5,0.25,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(1,6,0.25,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(5,10,0.25,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(50,250,10,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(8,15,0.5,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)
 