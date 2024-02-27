% Define X and Y
% 
% This script defines and assembles training (X) and target (Y, i.e., pCO2)
% variables for five US large Marine Ecosystems (Alaska, California
% Current, Insular Pacific / Hawaii, Gulf of Mexico / Caribbean, and US
% East Coast) in preparation for clustering via self-organizing maps and
% algorithm training.
% 
% Written by J.D. Sharp: 9/20/22
% Last updated by J.D. Sharp: 5/15/22
% 

%% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

for n = 1:length(region)

    %% load gridded pCO2 and predictors
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');

    %% Compute regional average total and seasonal data coverage
    disp(region{n});
    disp(['Total: ' num2str(mean(100.*SOCAT_grid.(region{n}).num_months(SOCAT_grid.(region{n}).idxspc(:,:,1))./SOCAT_grid.(region{n}).dim.z))]);
    disp(['Seasonal: ' num2str(mean(100.*SOCAT_grid.(region{n}).num_months_clim(SOCAT_grid.(region{n}).idxspc(:,:,1))./12))]);

    %% display status
    disp(['Processing data for clustering and algorithm training (' region{n} ')']);

    %% Create X and Y variables, as well as indices
    % define headers
    Vars_array.(region{n}).headers = ...
        {'Lon.' 'Lat.' 'month' 'sin_month' 'cos_month' 'year' 'SSS' 'SSH' ...
         'SST' 'IceC' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' 'MSLP' 'pCO2_atm'};
    Vars_array.(region{n}).headers_var = ...
        {'Lon.' 'Lat.' 'SSS' 'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'MLD' 'MSLP' 'pCO2_atm'};

    % create temporary versions of grid variables
    lon_tmp = repmat(Preds_grid.(region{n}).lon,1,Preds_grid.(region{n}).dim.y,...
        Preds_grid.(region{n}).dim.z);
    lat_tmp = repmat(Preds_grid.(region{n}).lat',Preds_grid.(region{n}).dim.x,1,...
        Preds_grid.(region{n}).dim.z);
    month_of_year_tmp = repmat(permute(Preds_grid.(region{n}).month_of_year,[3 2 1]),...
        Preds_grid.(region{n}).dim.x,Preds_grid.(region{n}).dim.y,1);
    sin_month_of_year_tmp = repmat(permute(sin((2.*pi.*Preds_grid.(region{n}).month_of_year)./12),[3 2 1]),...
        Preds_grid.(region{n}).dim.x,Preds_grid.(region{n}).dim.y,1);
    cos_month_of_year_tmp = repmat(permute(cos((2.*pi.*Preds_grid.(region{n}).month_of_year)./12),[3 2 1]),...
        Preds_grid.(region{n}).dim.x,Preds_grid.(region{n}).dim.y,1);
    year_tmp = repmat(permute(Preds_grid.(region{n}).year,[3 2 1]),...
        Preds_grid.(region{n}).dim.x,Preds_grid.(region{n}).dim.y);
    dist_tmp = repmat(Preds_grid.(region{n}).Dist,1,1,...
        Preds_grid.(region{n}).dim.z);
    bathy_tmp = repmat(Preds_grid.(region{n}).Bathy,1,1,...
        Preds_grid.(region{n}).dim.z);

    % create indices
    idx_clust = ~isnan(Preds_grid.(region{n}).SSS) & ...
           repmat(Preds_grid.(region{n}).idxspc(:,:,1),1,1,Preds_grid.(region{n}).dim.z);
    idx_clust_var = ~isnan(Preds_grid.(region{n}).SSS(:,:,1)) & ...
           Preds_grid.(region{n}).idxspc(:,:,1);
    idx_mod = ~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd) & ...
           repmat(Preds_grid.(region{n}).idxspc(:,:,1),1,1,Preds_grid.(region{n}).dim.z);

    % calculate variability
    Preds_grid.(region{n}).SSS_var = std(Preds_grid.(region{n}).SSS,[],3);
    Preds_grid.(region{n}).SSH_var = std(Preds_grid.(region{n}).SSH,[],3);
    Preds_grid.(region{n}).SST_var = std(Preds_grid.(region{n}).SST,[],3);
    Preds_grid.(region{n}).IceC_var = std(Preds_grid.(region{n}).IceC,[],3);
    Preds_grid.(region{n}).CHL_var = std(Preds_grid.(region{n}).CHL,[],3);
    Preds_grid.(region{n}).WindSpeed_var = std(Preds_grid.(region{n}).WindSpeed,[],3);
    Preds_grid.(region{n}).MLD_var = std(Preds_grid.(region{n}).MLD,[],3);
    Preds_grid.(region{n}).mslp_var = std(Preds_grid.(region{n}).mslp,[],3);
    Preds_grid.(region{n}).pCO2_atm_var = std(Preds_grid.(region{n}).pCO2_atm,[],3);
    lon_tmp_2 = mean(lon_tmp,3);
    lat_tmp_2 = mean(lat_tmp,3);

    % create X and Y
    Vars_array.(region{n}).X_clust = [lon_tmp(idx_clust) lat_tmp(idx_clust)...
         month_of_year_tmp(idx_clust) sin_month_of_year_tmp(idx_clust) ...
         cos_month_of_year_tmp(idx_clust) year_tmp(idx_clust) ...
         Preds_grid.(region{n}).SSS(idx_clust) Preds_grid.(region{n}).SSH(idx_clust) ...
         Preds_grid.(region{n}).SST(idx_clust) Preds_grid.(region{n}).IceC(idx_clust) ...
         log10(Preds_grid.(region{n}).CHL(idx_clust)) ...
         Preds_grid.(region{n}).WindSpeed(idx_clust) bathy_tmp(idx_clust) dist_tmp(idx_clust) ...
         Preds_grid.(region{n}).MLD(idx_clust) Preds_grid.(region{n}).mslp(idx_clust) ...
         Preds_grid.(region{n}).pCO2_atm(idx_clust)];
    Vars_array.(region{n}).X_clust_var = [
         lon_tmp_2(idx_clust_var) lat_tmp_2(idx_clust_var) ...
         Preds_grid.(region{n}).SSS_var(idx_clust_var) Preds_grid.(region{n}).SSH_var(idx_clust_var) ...
         Preds_grid.(region{n}).SST_var(idx_clust_var) Preds_grid.(region{n}).IceC_var(idx_clust_var) ...
         log10(Preds_grid.(region{n}).CHL_var(idx_clust_var)) ...
         Preds_grid.(region{n}).WindSpeed_var(idx_clust_var) ...
         Preds_grid.(region{n}).MLD_var(idx_clust_var) Preds_grid.(region{n}).mslp_var(idx_clust_var) ...
         Preds_grid.(region{n}).pCO2_atm_var(idx_clust_var)];
    Vars_array.(region{n}).X_mod = [lon_tmp(idx_mod) lat_tmp(idx_mod) month_of_year_tmp(idx_mod) sin_month_of_year_tmp(idx_mod) ...
         cos_month_of_year_tmp(idx_mod) year_tmp(idx_mod) ...
         Preds_grid.(region{n}).SSS(idx_mod) Preds_grid.(region{n}).SSH(idx_mod) ...
         Preds_grid.(region{n}).SST(idx_mod) Preds_grid.(region{n}).IceC(idx_mod)...
         log10(Preds_grid.(region{n}).CHL(idx_mod)) ...
         Preds_grid.(region{n}).WindSpeed(idx_mod) bathy_tmp(idx_mod) dist_tmp(idx_mod) ...
         Preds_grid.(region{n}).MLD(idx_mod) Preds_grid.(region{n}).mslp(idx_mod) ...
         Preds_grid.(region{n}).pCO2_atm(idx_mod)];
    Vars_array.(region{n}).Y_mod = SOCAT_grid.(region{n}).fco2_ave_wtd(idx_mod);

    % add indices to gridded predictors
    Preds_grid.(region{n}).idx_clust = idx_clust;
    Preds_grid.(region{n}).idx_clust_var = idx_clust_var;
    Preds_grid.(region{n}).idx_mod = idx_mod;

    % Save predictor/target arrays and gridded predictor data with cluster indices
    save(['Data/' region{n} '/variable_arrays'],'Vars_array','-v7.3');
    save(['Data/' region{n} '/gridded_predictors'],'Preds_grid','-v7.3');

    % clean up
    clear lon_tmp lat_tmp month_of_year_tmp cos_month_of_year_tmp sin_month_of_year_tmp
    clear year_tmp dist_tmp bathy_tmp idx_clust idx_clust_var idx_mod
    clear SOCAT_grid Vars_array Preds_grid

end

clear n
