function create_figures(vrs,lme_shape,lme_idx,region)

for n = 1:length(region)

    % load RFR-LME
    load(['Data/RFR-LME/' vrs '_' region{n}],'RFR_LME');

    % SSS
    plot_temporal_mean(RFR_LME.lim,RFR_LME.dim,RFR_LME.lat,RFR_LME.lon,RFR_LME.SSS,...
          cmocean('haline',20),'SSS','Sea Surface Salinity',region{n},lme_shape,lme_idx);
    % SST
    plot_temporal_mean(RFR_LME.lim,RFR_LME.dim,RFR_LME.lat,RFR_LME.lon,RFR_LME.SST,...
          cmocean('thermal',20),'SST','Sea Surface Temperature',region{n},lme_shape,lme_idx);

end
