% this script defines the bounds of the eleven LMEs
define_regions_eiwg
% plot surface predictor variables across full region
plot_temporal_mean_full(29.75,38.25,cmocean('haline'),'SSS','Sea Surface Salinity',region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.025,0.105,parula,'SSH','Sea Surface Height Anomaly (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(1,31,cmocean('thermal'),'SST',['Sea Surface Temperature (' char(176) 'C)'],region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.025,1.025,cmocean('tempo'),'IceC','Sea Ice Concentration Fraction',region,lme_shape,lme_idx)
plot_temporal_mean_full(-0.025,1.025,cmocean('algae'),'CHL','Sea Surface Chlorophyll (log_{10})',region,lme_shape,lme_idx)
plot_temporal_mean_full(4.5,10.5,cmocean('amp'),'WindSpeed','Wind Speed (m s^{-1})',region,lme_shape,lme_idx)
plot_temporal_mean_full(-125,6125,cmocean('deep'),'Bathy','Bottom Depth (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(-2.5,52.5,jet,'MLD','Mixed Layer Depth (m)',region,lme_shape,lme_idx)
plot_temporal_mean_full(0.990,1.010,cmocean('dense'),'mslp','Sea Level Pressure (atm)',region,lme_shape,lme_idx)
plot_temporal_mean_full(369.5,390.5,cmocean('solar'),'pCO2_atm','Atmospheric pCO_{2} (\muatm)',region,lme_shape,lme_idx)
