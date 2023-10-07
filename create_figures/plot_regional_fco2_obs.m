% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

% extract each LME from large grid
for n = 1:length(region)
    % load data
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    % plot detrended gridded mean pCO2
    plot_temporal_mean(SOCAT_grid.(region{n}).lim,...
        SOCAT_grid.(region{n}).dim,SOCAT_grid.(region{n}).lat,...
        SOCAT_grid.(region{n}).lon,SOCAT_grid.(region{n}).fco2_ave_wtd_detrend,...
        parula(20),'fCO2_obs','Surface {\itf}CO_{2} (\muatm)',region{n});
end
