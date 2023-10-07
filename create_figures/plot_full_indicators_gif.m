%% create region-wide animations of RFR-LME indicators

plot_full_gif(295,475,parula,'fCO2','Sea Surface {\itf}_{CO2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx);
plot_full_gif(295,475,parula,'pCO2','Sea Surface {\itp}CO_{2(RFR-LME)} (\muatm)',region,lme_shape,lme_idx);

        plot_full_gif(2050,2400,flipud(cmocean('deep')),'TA','Sea Surface {\itA}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
        plot_full_gif(1850,2200,cmocean('matter'),'DIC','Sea Surface {\itC}_{T(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
        plot_full_gif(8.0,8.2,cmocean('solar'),'pH','Sea Surface pH_{T(RFR-LME)}',region,lme_shape,lme_idx)
        plot_full_gif(0.5,4.5,flipud(cmocean('ice')),'OmA','Sea Surface \Omega_{A(RFR-LME)}',region,lme_shape,lme_idx)
        plot_full_gif(1.0,6.0,flipud(cmocean('ice')),'OmC','Sea Surface \Omega_{C(RFR-LME)}',region,lme_shape,lme_idx)
        plot_full_gif(5,10,flipud(cmocean('solar')),'H','Sea Surface [H^{+}]_{T(RFR-LME)} (nmol kg^{-1})',region,lme_shape,lme_idx)
        plot_full_gif(50,250,flipud(cmocean('ice')),'CO3','Sea Surface [CO_{3}^{2-}]_{(RFR-LME)} (\mumol kg^{-1})',region,lme_shape,lme_idx)
        plot_full_gif(9,16,cmocean('speed'),'RF','Sea Surface RF_{(RFR-LME)}',region,lme_shape,lme_idx)
        plot_full_gif(0,40,cmocean('amp',21),'ufCO2','Sea Surface {\itf}CO_{2(RFR-LME)} Uncertainty',region,lme_shape,lme_idx)
   