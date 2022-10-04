% Fit Algorithms
% 
% This script trains neural networks and random forest regressions for the
% prediction of sea surface fCO2 in various clusters for five US large
% Marine Ecosystems (Alaska, California Current, Insular Pacific / Hawaii,
% Gulf of Mexico / Caribbean, and US East Coast).
% 
% Written by J.D. Sharp: 8/26/22
% Last updated by J.D. Sharp: 9/28/22
% 

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 1:length(region)

    %% load gridded pCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');

    %% define variables to use for models
    Vars_array.(region{n}).vars = ...
        {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
            'SSH' 'SST' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
            'MSLP' 'pCO2_atm'};
    Vars_array.(region{n}).idx_vars = ...
        find(ismember(Vars_array.(region{n}).headers,...
        Vars_array.(region{n}).vars));

    %% define group indices
    Vars_array.(region{n}).idx_group_mod = ...
        Preds_grid.(region{n}).groups(Preds_grid.(region{n}).idx_mod);

    %% fit NNs
    disp(['Fitting NNs (' region{n} ')']);
    % fit model
    size1 = [30 25 20];
    size2 = [10 15 20];
    [Mods.(region{n}).nn,Mods.(region{n}).err_nn,...
     Mods.(region{n}).rmse_nn,Mods.(region{n}).r2_nn,...
     Mods.(region{n}).Y_fit_nn,Mods.(region{n}).delta_nn] = ...
        fit_nn(Vars_array.(region{n}).X_mod,Vars_array.(region{n}).Y_mod,...
        Vars_array.(region{n}).idx_vars,Vars_array.(region{n}).idx_group_mod,size1,size2);
    
    % plot delta pco2 for NNs
    x_edges = 0:5:1000; x_mids = 2.5:5:997.5;
    y_edges = -1000:10:1000; y_mids = -995:10:995;
    plot_delta(x_edges,x_mids,y_edges,y_mids,Mods.(region{n}).Y_fit_nn.all,...
        Mods.(region{n}).delta_nn.all,region{n},[200 800],[-300 300],'NN');
    clear x_edges x_mids y_edges y_mids

    %% fit RFRs
    disp(['Fitting RFRs (' region{n} ')']);
    % set parameters for random forest
    nTrees = 1000;
    minLeafSize = 5;
    numpredictors = 6;
    % fit model
    [Mods.(region{n}).rfr,Mods.(region{n}).err_rfr,...
     Mods.(region{n}).rmse_rfr,Mods.(region{n}).r2_rfr,...
     Mods.(region{n}).Y_fit_rfr,Mods.(region{n}).delta_rfr] = ...
        fit_rfr(Vars_array.(region{n}).X_mod,Vars_array.(region{n}).Y_mod,...
        Vars_array.(region{n}).idx_vars,Vars_array.(region{n}).idx_group_mod,...
        Vars_array.(region{n}).headers,nTrees,minLeafSize,numpredictors);
    
    % plot delta pco2 for RFRs
    x_edges = 0:5:1000; x_mids = 2.5:5:997.5;
    y_edges = -1000:10:1000; y_mids = -995:10:995;
    plot_delta(x_edges,x_mids,y_edges,y_mids,Mods.(region{n}).Y_fit_rfr.all,...
        Mods.(region{n}).delta_rfr.all,region{n},[200 800],[-300 300],'RFR');
    clear x_edges x_mids y_edges y_mids

    % clean up
    clear nTrees minLeafSize numpredictors

    % Save models and predictor/target arrays with variable indices
    save(['Data/' region{n} '/us_lme_models'],'Mods','-v7.3');
    save(['Data/' region{n} '/variable_arrays'],'Vars_array','-v7.3');
    clearvars -except region

end
