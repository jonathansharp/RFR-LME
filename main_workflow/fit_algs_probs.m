% Fit Algorithms
% 
% This script trains neural networks and random forest regressions for the
% prediction of sea surface fCO2 in various US LMEs.
% 
% Written by J.D. Sharp: 8/26/22
% Last updated by J.D. Sharp: 2/16/23
% 

for n = 1:length(region)

    %% load gridded pCO2, predictors, and clusters
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    if test_idx == 0
        load(['Data/' region{n} '/gridded_clusters'],'Clusts_grid');
    elseif test_idx == 1
        load(['Data/' region{n} '/gridded_clusters_test'],'Clusts_grid');
    end

    %% define variables to use for models
    % exclude Chlorophyll and SSH as predictors for Bering-Chukchi and Beaufort Seas
    if n == 5 || n == 6
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SST' 'IceC' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
    elseif n == 2 || n == 3 || n == 4
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
    % exclude Sea Ice Concentration as predictors for south of AK
    else
        Vars_array.(region{n}).vars = ...
            {'Lon.' 'Lat.' 'sin_month' 'cos_month' 'year' 'SSS' ...
                'SSH' 'SST' 'CHL' 'Wind' 'Bathy.' 'Dist.' 'MLD' ...
                'MSLP' 'pCO2_atm'};
    end
    % determine variable index
    Vars_array.(region{n}).idx_vars = ...
        find(ismember(Vars_array.(region{n}).headers,...
        Vars_array.(region{n}).vars));

    %% define group indices
    % pre-allocate
    Vars_array.(region{n}).prob_group_mod = ...
        nan(length(Vars_array.(region{n}).Y_mod),num_groups(n));
    Vars_array.(region{n}).idx_group_mod = ...
        Clusts_grid.(region{n}).groups(Preds_grid.(region{n}).idx_mod);
    for c = 1:num_groups(n)
        Vars_array.(region{n}).prob_group_mod(:,c) = ...
            Clusts_grid.(region{n}).probabilities.(['c' num2str(c)])(Preds_grid.(region{n}).idx_mod) > 0.1;
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
        Vars_array.(region{n}).idx_vars,Vars_array.(region{n}).idx_group_mod,...
        Vars_array.(region{n}).prob_group_mod,Vars_array.(region{n}).headers,nTrees,minLeafSize,...
        numpredictors,num_groups(n),region{n});
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
        Val.(region{n}).delta_rfr.all;
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
    plot_delta(x_edges,x_mids,y_edges,y_mids,Val.(region{n}).Y_fit_rfr.all,...
        Val.(region{n}).delta_rfr.all,region{n},[200 800],[-300 300],'RFR');

    % clean up
    clear x_edges x_mids y_edges y_mids

    %% Plot gridded delta pCO2 for RFRs
    figure('visible','off');
    worldmap([Preds_grid.(region{n}).lim.latmin ...
        Preds_grid.(region{n}).lim.latmax],...
       [Preds_grid.(region{n}).lim.lonmin ...
        Preds_grid.(region{n}).lim.lonmax]);
    set(gca,'fontsize',16);
    pcolorm(repmat(Preds_grid.(region{n}).lat',Preds_grid.(region{n}).dim.x,1),...
            repmat(Preds_grid.(region{n}).lon,1,Preds_grid.(region{n}).dim.y),...
            mean(abs(Val.(region{n}).delta_rfr_grid),3,'omitnan'));
    plot_land('map');
    c=colorbar;
    colormap(cmocean('tempo',16));
    caxis([-1 31]);
    c.TickLength = 0;
    c.Label.String = '\Delta{\itf}CO_{2} (\muatm)';
    cbarrow('up');
    clear c

    % save figure
    if ~isfolder('Figures'); mkdir('Figures'); end
    exportgraphics(gcf,['Figures/' region{n} '_delta_fCO2_RFR_gridded.png']);
    close

    %% Save models and predictor/target arrays with variable indices
    if test_idx == 0
        save(['Data/' region{n} '/us_lme_models'],'Mods','-v7.3');
        save(['Data/' region{n} '/us_lme_model_evals'],'Val','-v7.3');
    elseif test_idx == 1
        save(['Data/' region{n} '/us_lme_models_test'],'Mods','-v7.3');
        save(['Data/' region{n} '/us_lme_model_evals_test'],'Val','-v7.3');
    end
    save(['Data/' region{n} '/variable_arrays'],'Vars_array','-v7.3');

    %% clean up
    clear Mods Vars_array

end

%% plot gridded delta values across region
cmap_type = 'cmocean';
cmap_name = 'balance';
zero_piv = 1;
cmap_segs = 9;
plot_delta_mean_full(-22.5,22.5,cmap_type,cmap_name,cmap_segs,zero_piv,...
    num_groups,'delta_rfr_grid','\Delta{\itf}CO_{2}',region,lme_shape,lme_idx);

%% plot gridded absolute delta values across region
cmap_type = 'cmocean';
cmap_name = 'amp';
zero_piv = 0;
cmap_segs = 9;
plot_delta_mean_full(-1.25,21.25,cmap_type,cmap_name,cmap_segs,...
    zero_piv,num_groups,'delta_rfr_grid_abs','\Delta{\itf}CO_{2}',...
    region,lme_shape,lme_idx);

%% plot delta values for each LME
% fCO2_rfr_errors
