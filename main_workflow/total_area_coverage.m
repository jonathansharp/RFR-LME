% total area coverage
define_regions_eiwg

for n = 1:length(region)

    %% load gridded pCO2 and predictors
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');

    %% Compute regional average total and seasonal data coverage
    disp(region{n});
    total_coverage(n) = ...
        mean(100.*SOCAT_grid.(region{n}).num_months(SOCAT_grid.(region{n}).idxspc(:,:,1))./SOCAT_grid.(region{n}).dim.z);
    seasonal_coverage(n) = mean(100.*SOCAT_grid.(region{n}).num_months_clim(SOCAT_grid.(region{n}).idxspc(:,:,1))./12);
    area(n) = SOCAT_grid.(region{n}).area_tot_km2;

end