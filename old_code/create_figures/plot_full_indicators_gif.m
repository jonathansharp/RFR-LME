%% create region-wide animations of RFR-LME indicators
addpath(genpath(pwd));
% this script defines the bounds of the eleven LMEs
define_regions_eiwg

plot_full_gif(295,475,flipud(slanCM('romao')),'fco2_ave_wtd','Sea Surface {\itf}_{CO2(SOCAT)} (\muatm)',region,lme_shape,lme_idx);
plot_full_gif(295,475,flipud(slanCM('romao')),'fCO2','Sea Surface {\itf}_{CO2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx);
plot_full_gif(295,475,flipud(slanCM('romao')),'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx);
plot_full_gif(2000,2400,flipud(slanCM('romao')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_full_gif(1800,2100,flipud(slanCM('romao')),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_full_gif(8.0,8.2,flipud(slanCM('romao')),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
plot_full_gif(0,5.0,flipud(slanCM('romao')),'OmA','Sea Surface \Omega_{A(RFR-LME)}',region,lme_shape,lme_idx)
plot_full_gif(0,7.0,flipud(slanCM('romao')),'OmC','Sea Surface \Omega_{C(RFR-LME)}',region,lme_shape,lme_idx)
plot_full_gif(7,10,flipud(slanCM('romao')),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
plot_full_gif(0,300,flipud(slanCM('romao')),'CO3','Sea Surface [CO_{3}^{2-}]_{(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
plot_full_gif(6,18,flipud(slanCM('romao')),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)

%plot_full_gif(0,40,cmocean('amp',21),'ufCO2','Sea Surface {\itf}CO_{2(RFR-LME)} Uncertainty',region,lme_shape,lme_idx)
