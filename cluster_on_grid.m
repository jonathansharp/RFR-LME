% Cluster data
% 
% This script uses the specified satellite and reanalysis variables for
% five US large Marine Ecosystems (Alaska, California Current, Insular
% Pacific / Hawaii, Gulf of Mexico / Caribbean, and US East Coast) to
% define time-varying clusters within which to train machine learning
% algorithms via a self-organizing map method.
% 
% Written by J.D. Sharp: 8/26/22
% Last updated by J.D. Sharp: 9/20/22
% 

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 1:length(region)

    %% load gridded pCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');

    %% cluster gridded data
    disp(['Clustering data (' region{n} ')']);
    % define variable index
    vars = {'SSS' 'SSH' 'SST' 'Wind'};
    idx_vars = find(ismember(Vars_array.(region{n}).headers,vars));
    % obtain gridded group numbers
    Preds_grid.(region{n}).groups = ...
        cluster_data(Vars_array.(region{n}).X_clust,idx_vars,...
        Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
        Preds_grid.(region{n}).month_of_year,...
        Preds_grid.(region{n}).idx_clust,region{n});
    % clean up
    clear vars idx_vars idx_clust

    % Save gridded predictor data (with clusters included) for all LMEs
    save(['Data/' region{n} '/gridded_predictors'],'Preds_grid','-v7.3');
    clearvars -except region

end
