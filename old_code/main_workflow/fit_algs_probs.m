% Fit Algorithms
% 
% This script trains random forest regressions for the prediction of sea
% surface fCO2 in various US LMEs.
% 
% Written by J.D. Sharp: 8/26/22
% Last updated by J.D. Sharp: 10/17/24
% 

% include moorings?
if include_moorings == 0; exts = '_no_moorings'; else; exts = ''; end

for n = 1:length(region)

    %% load gridded pCO2, predictors, and clusters
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays' exts],'Vars_array');
    if gmm_test_idx == 0
        load(['Data/' region{n} '/gridded_clusters'],'Clusts_grid');
    elseif gmm_test_idx == 1
        load(['Data/' region{n} '/gridded_clusters_test'],'Clusts_grid');
    end

    %% define variables to use for models
    % exclude Chlorophyll and SSH as predictors for Bering-Chukchi and Beaufort Seas
    if n == 5 || n == 6
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SST' 'IceC' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
%         Vars_array.(region{n}).vars = ...
%             {'SSS' ...
%                 'SST' 'IceC' 'Wind' 'Bathy.' 'MLD' ...
%                 'MSLP' 'pCO2_atm'};
    elseif n == 2 || n == 3 || n == 4 % Aleutian, GoA, EBS
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
%         Vars_array.(region{n}).vars = ...
%             {'SSS' ...
%                 'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'Bathy.' 'MLD' ...
%                 'MSLP' 'pCO2_atm'};
    % exclude Sea Ice Concentration as predictors for south of AK
    else
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SSH' 'SST' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
%         Vars_array.(region{n}).vars = ...
%             {'SSS' ...
%                 'SSH' 'SST' 'CHL' 'Wind' 'Bathy.' 'MLD' ...
%                 'MSLP' 'pCO2_atm'};
    end
    % determine variable index
    Vars_array.(region{n}).idx_vars = ...
        find(ismember(Vars_array.(region{n}).headers,...
        Vars_array.(region{n}).vars));

    %% define group indices
    % pre-allocate
    Vars_array.(region{n}).idx_group_mod = ...
        Clusts_grid.(region{n}).groups(Preds_grid.(region{n}).idx_mod);
    Vars_array.(region{n}).prob_group_mod = ...
        nan(length(Vars_array.(region{n}).Y_mod),num_groups(n));
    Vars_array.(region{n}).probs_groups = ...
        nan(length(Vars_array.(region{n}).Y_mod),num_groups(n));
    for c = 1:num_groups(n)
        Vars_array.(region{n}).prob_group_mod(:,c) = ...
            Clusts_grid.(region{n}).probabilities.(['c' num2str(c)])(Preds_grid.(region{n}).idx_mod) > 0.10;
        Vars_array.(region{n}).probs_groups(:,c) = ...
            Clusts_grid.(region{n}).probabilities.(['c' num2str(c)])(Preds_grid.(region{n}).idx_mod);
    end

%     %% fit NNs
%     rng(7); % for reproducibility
%     disp(['Fitting NNs (' region{n} ')']);
%     % fit model
%     size1 = [30 25 20];
%     size2 = [10 15 20];
%     [Mods.(region{n}).nn,Mods.(region{n}).err_nn,...
%      Mods.(region{n}).rmse_nn,Mods.(region{n}).r2_nn,...
%      Mods.(region{n}).Y_fit_nn,Mods.(region{n}).delta_nn] = ...
%         fit_nn(Vars_array.(region{n}).X_mod,Vars_array.(region{n}).Y_mod,...
%         Vars_array.(region{n}).idx_vars,Vars_array.(region{n}).idx_group_mod,size1,size2);
%     
%     %% plot 2D histogram of delta pco2 for NNs
%     x_edges = 0:5:1000; x_mids = 2.5:5:997.5;
%     y_edges = -1000:10:1000; y_mids = -995:10:995;
%     plot_delta(x_edges,x_mids,y_edges,y_mids,Mods.(region{n}).Y_fit_nn.all,...
%         Mods.(region{n}).delta_nn.all,region{n},[200 800],[-300 300],'NN');
%     clear x_edges x_mids y_edges y_mids

    %% fit RFRs
    rng(5); % for reproducibility
%     minLeafSize = 5:5:50;
%     test_rmse = nan(size(minLeafSize));
%     colors = 'rbcmygrbcm';
%     figure;
%     hold on;
%     for lf = 1:length(minLeafSize)
    disp(['Fitting RFRs (' region{n} ')']);
    % fit model
    [Mods.(region{n}).rfr,Val.(region{n}).avg_err_rfr,Val.(region{n}).med_err_rfr,...
     Val.(region{n}).rmse_rfr,Val.(region{n}).mae_rfr,Val.(region{n}).r2_rfr,...
     Val.(region{n}).Y_fit_rfr,Val.(region{n}).delta_rfr] = ...
        fit_rfr(Vars_array.(region{n}).X_mod,Vars_array.(region{n}).Y_mod,...
        Vars_array.(region{n}).idx_vars,Vars_array.(region{n}).probs_groups,...
        Vars_array.(region{n}).prob_group_mod,Vars_array.(region{n}).headers,...
        nTrees,minLeafSize,ceil(length(Vars_array.(region{n}).idx_vars)/3),...
        num_groups(n),region{n});
    % print error stats
    disp(['Avg. Err. = ' num2str(round(Val.(region{n}).avg_err_rfr(end),2))]);
    disp(['Med. Err. = ' num2str(round(Val.(region{n}).med_err_rfr(end),2))]);
    disp(['RMSE = ' num2str(round(Val.(region{n}).rmse_rfr(end),2))]);
    disp(['R2 = ' num2str(round(Val.(region{n}).r2_rfr(end),2))]);
%     test_rmse(lf) = Mods.(region{n}).rmse_rfr(2);
    % place errors on grid
    Val.(region{n}).delta_rfr_grid = ...
        nan(size(Preds_grid.(region{n}).idx_mod));
    Val.(region{n}).delta_rfr_grid(Preds_grid.(region{n}).idx_mod) = ...
        Val.(region{n}).delta_rfr.all(:,end);
    Val.(region{n}).delta_rfr_grid_abs = ...
        abs(Val.(region{n}).delta_rfr_grid);
    % add dimensions to grid
    Val.(region{n}).lat = Preds_grid.(region{n}).lat;
    Val.(region{n}).lon = Preds_grid.(region{n}).lon;
    Val.(region{n}).month = Preds_grid.(region{n}).month;
%     end
%     xlabel('Number of Grown Trees');
%     ylabel('Root Mean Squared Error');
%     legend({'5' '10' '15' '20' '25' '30' '35' '40' '45' '50'});
%     exportgraphics(gcf,['Figures/' region{n} '/LeafSizeTest_zoom.png']);

    %% plot 2D histogram of delta pco2 for RFRs
    x_edges = 0:5:1000;
    x_mids = 2.5:5:997.5;
    y_edges = -1000:10:1000;
    y_mids = -995:10:995;
    plot_delta(x_edges,x_mids,y_edges,y_mids,Val.(region{n}).Y_fit_rfr.all(:,end),...
        Val.(region{n}).delta_rfr.all(:,end),region{n},[200 800],[-300 300],'RFR');

    % clean up
    clear x_edges x_mids y_edges y_mids

    %% plot 2D histogram of pco2 vs. pco2 for RFRs
    x_edges = 0:5:1000;
    x_mids = 2.5:5:997.5;
    y_edges = 0:5:1000;
    y_mids = 0:5:1000;
    plot_mod_v_meas(x_edges,x_mids,y_edges,y_mids,Vars_array.(region{n}).Y_mod,...
        Val.(region{n}).Y_fit_rfr.all(:,end), Val.(region{n}).delta_rfr.all(:,end),...
        region{n},[200 600],[200 600],'RFR');

    % clean up
    clear x_edges x_mids y_edges y_mids

    %% Save models and predictor/target arrays with variable indices
    if rfr_test_idx == 0
        save(['Data/' region{n} '/us_lme_models' exts],'Mods','-v7.3');
        save(['Data/' region{n} '/us_lme_model_evals' exts],'Val','-v7.3');
    elseif rfr_test_idx == 1
        save(['Data/' region{n} '/us_lme_models_test' exts],'Mods','-v7.3');
        save(['Data/' region{n} '/us_lme_model_evals_test' exts],'Val','-v7.3');
    end
    save(['Data/' region{n} '/variable_arrays' exts],'Vars_array','-v7.3');

    %% clean up
    clear Mods Val Vars_array

end
