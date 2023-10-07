% create region-wide long-term mean plots of RFR-LME indicators


% plot fCO2 and pCO2 across full region
plot_temporal_mean_full(295,475,parula,'fCO2','Sea Surface {\itf}_{CO2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(295,475,parula,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
 
% plot OA indicators trends across full region
plot_temporal_trend_full(-5,5,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} Trend (\muatm yr^{-1})',region,lme_shape,lme_idx)

% plot fCO2 and pCO2 across full region (seasonally)
plot_temporal_mean_full_seas(295,475,parula,'fCO2','Sea Surface {\itf}_{CO2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)
plot_temporal_mean_full_seas(295,475,parula,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx)


    %% plot OA indicators across full region
    plot_temporal_mean_full(2050,2400,flipud(cmocean('deep')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full(1850,2200,cmocean('matter'),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full(8.0,8.2,cmocean('solar'),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full(0.5,5.5,flipud(cmocean('algae')),'OmA','Sea Surface \Omega_{ar(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full(1.0,7.0,flipud(cmocean('dense')),'OmC','Sea Surface \Omega_{ca(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full(5,9.5,cmocean('ice'),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full(50,300,cmocean('turbid'),'CO3','Sea Surface [CO_{3}^{2-}]_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full(8,16,cmocean('speed'),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full(0,50,cmocean('amp'),'upCO2','Sea Surface {\itp}CO_{2(RFR-LME)} Uncertainty (\muatm)',region,lme_shape,lme_idx)
    plot_temporal_mean_full(0,50,cmocean('amp'),'ufCO2','Sea Surface {\itf}_{CO2(RFR-LME)} Uncertainty (\muatm)',region,lme_shape,lme_idx)
    plot_temporal_mean_full(1,1.2,flipud(cmocean('tempo')),'TA_DIC','Sea Surface TA/DIC_{(RFR-LME)}',region,lme_shape,lme_idx)
    
    %% plot OA indicators trends across full region
    plot_temporal_trend_full(2050,2400,flipud(cmocean('deep')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    
        %% plot OA indicators across full region (seasonally)
    plot_temporal_mean_full_seas(2050,2400,flipud(cmocean('deep')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(1850,2200,cmocean('matter'),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(8.0,8.2,cmocean('solar'),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(0.5,4.5,flipud(cmocean('ice')),'OmA','Sea Surface \Omega_{ar(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(1.0,6.0,flipud(cmocean('ice')),'OmC','Sea Surface \Omega_{ca(RFR-LME)}',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(5,10,flipud(cmocean('solar')),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(50,250,flipud(cmocean('ice')),'CO3','Sea Surface [CO_{3}^{2-}]_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
    plot_temporal_mean_full_seas(9,16,cmocean('speed'),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)
