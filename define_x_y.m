% Define X and Y
% 
% This script defines and assembles training (X) and target (Y, i.e., pCO2)
% variables for five US large Marine Ecosystems (Alaska, California
% Current, Insular Pacific / Hawaii, Gulf of Mexico / Caribbean, and US
% East Coast) in preparation for clustering via self-organizing maps and
% algorithm training.
% 
% Written by J.D. Sharp: 9/20/22
% Last updated by J.D. Sharp: 9/20/22
% 

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 1:length(region)

    %% load gridded pCO2 and predictors
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');

    %% display status
    disp(['Processing data for clustering and algorithm training (' region{n} ')']);

    %% Create X and Y variables, as well as indices
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
           repmat(Preds_grid.(region{n}).ocean_mask,1,1,Preds_grid.(region{n}).dim.z);
    idx_mod = ~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd) & ...
           repmat(Preds_grid.(region{n}).ocean_mask,1,1,Preds_grid.(region{n}).dim.z);

    % create X and Y
    Vars_array.(region{n}).X_clust = [lon_tmp(idx_clust) lat_tmp(idx_clust) month_of_year_tmp(idx_clust) sin_month_of_year_tmp(idx_clust) ...
         cos_month_of_year_tmp(idx_clust) year_tmp(idx_clust) ...
         Preds_grid.(region{n}).SSS(idx_clust) Preds_grid.(region{n}).SSH(idx_clust) ...
         Preds_grid.(region{n}).SST(idx_clust) log10(Preds_grid.(region{n}).CHL(idx_clust)) ...
         Preds_grid.(region{n}).WindSpeed(idx_clust) bathy_tmp(idx_clust) dist_tmp(idx_clust) ...
         Preds_grid.(region{n}).MLD(idx_clust) Preds_grid.(region{n}).mslp(idx_clust) ...
         Preds_grid.(region{n}).pCO2_atm(idx_clust)];
    Vars_array.(region{n}).X_mod = [lon_tmp(idx_mod) lat_tmp(idx_mod) month_of_year_tmp(idx_mod) sin_month_of_year_tmp(idx_mod) ...
         cos_month_of_year_tmp(idx_mod) year_tmp(idx_mod) ...
         Preds_grid.(region{n}).SSS(idx_mod) Preds_grid.(region{n}).SSH(idx_mod) ...
         Preds_grid.(region{n}).SST(idx_mod) log10(Preds_grid.(region{n}).CHL(idx_mod)) ...
         Preds_grid.(region{n}).WindSpeed(idx_mod) bathy_tmp(idx_mod) dist_tmp(idx_mod) ...
         Preds_grid.(region{n}).MLD(idx_mod) Preds_grid.(region{n}).mslp(idx_mod) ...
         Preds_grid.(region{n}).pCO2_atm(idx_mod)];
    Vars_array.(region{n}).Y_mod = SOCAT_grid.(region{n}).fco2_ave_wtd(idx_mod);

    % save indices
    Preds_grid.(region{n}).idx_clust = idx_clust;
    Preds_grid.(region{n}).idx_mod = idx_mod;

    % headers for training data
    Vars_array.(region{n}).headers = ...
        {'Lon.' 'Lat.' 'month' 'sin_month' 'cos_month' 'year' 'SSS' ...
         'SSH' 'SST' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' 'MSLP' 'pCO2_atm'};

    % clean up
    clear lon_tmp lat_tmp cos_month_of_year_tmp sin_month_of_year_tmp
    clear month_of_year_tmp year_tmp dist_tmp bathy_tmp idx_clust idx_mod

    % Save predictor/target arrays and gridded predictor data with cluster indices
    save(['Data/' region{n} '/variable_arrays'],'Vars_array','-v7.3');
    save(['Data/' region{n} '/gridded_predictors'],'Preds_grid','-v7.3');
    clearvars -except region

end
