% create region-wide long-term mean plots of RFR-LME indicators

addpath(genpath(pwd));

% this script defines the bounds of the eleven LMEs
define_regions_eiwg

% plot fCO2 and pCO2 across full region
plot_temporal_mean_full(295,475,flipud(slanCM('romao')),'fCO2','Sea Surface {\itf}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(295,475,flipud(slanCM('romao')),'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
 
% plot pCO2 trends across full region
plot_temporal_trend_full(-5,5,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} Trend (\muatm yr^{-1})',region,lme_shape,lme_idx)

% plot pCO2 IAV across full region
plot_temporal_iav_full(0,30,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} IAV (\muatm)',region,lme_shape,lme_idx)

% plot fCO2 and pCO2 across full region (seasonally)
plot_temporal_mean_full_seas(295,475,parula,'fCO2','Sea Surface {\itf}_{CO2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(250,500,parula,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)


%% plot OA indicators across full region
plot_temporal_mean_full(2000,2400,flipud(slanCM('romao')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(1800,2100,flipud(slanCM('romao')),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(8.0,8.2,flipud(slanCM('romao')),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,5.0,flipud(slanCM('romao')),'OmA','Sea Surface \Omega_{ar(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,7.0,flipud(slanCM('romao')),'OmC','Sea Surface \Omega_{ca(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full(7,10,flipud(slanCM('romao')),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,300,flipud(slanCM('romao')),'CO3','Sea Surface [CO_{3}^{2-}]_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(6,18,flipud(slanCM('romao')),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,50,cmocean('amp'),'upCO2','Sea Surface {\itp}CO_{2(RFR-LME)} Uncer. (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,50,cmocean('amp'),'ufCO2','Sea Surface {\itf}CO_{2(RFR-LME)} Uncer. (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,0.1,cmocean('amp'),'uRF','Sea Surface RF Uncer.',region,lme_shape,lme_idx)
plot_temporal_mean_full(1.05,1.25,flipud(slanCM('romao')),'TA_DIC','Sea Surface TA/DIC_{(RFR-LME)}',region,lme_shape,lme_idx)

%% plot OA indicators trends across full region
%plot_temporal_trend_full(2050,2400,flipud(cmocean('deep')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)

%% plot OA indicators across full region (seasonally)
plot_temporal_mean_full_seas(2000,2400,flipud(slanCM('romao')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(1800,2100,flipud(slanCM('romao')),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(8.0,8.2,flipud(slanCM('romao')),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(0,5.0,flipud(slanCM('romao')),'OmA','Sea Surface \Omega_{ar(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(0,7.0,flipud(slanCM('romao')),'OmC','Sea Surface \Omega_{ca(RFR-LME)}',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(7,10,flipud(slanCM('romao')),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(50,300,flipud(slanCM('romao')),'CO3','Sea Surface [CO_{3}^{2-}]_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(6,18,flipud(slanCM('romao')),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)
