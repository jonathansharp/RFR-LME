% Predict fCO2 on grid
% 
% This script uses neural networks and random forest regressions to predict
% fCO2 in time-varying clusters defined by a self-organizing map method for
% US large Marine Ecosystems.
% 
% Written by J.D. Sharp: 9/15/22
% Last updated by J.D. Sharp: 2/23/22
% 

for n = 1:length(region)

    %% display status
    disp(['Predicting fCO2 (' region{n} ')']);

    %% load gridded fCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    load(['Data/' region{n} '/gridded_clusters'],'Clusts_grid');
    load(['Data/' region{n} '/us_lme_models'],'Mods');

    %% pre-allocate results
    fco2_rfr_tmp = nan(length(Vars_array.(region{n}).X_clust),num_groups(n));
    fco2_gmm_probs = nan(length(Vars_array.(region{n}).X_clust),num_groups(n));

    %% predict fCO2 using RFR for each cluster
    for c = 1:num_groups(n)

        % apply random forest regression to data
        temp_rfr = Mods.(region{n}).rfr.(['c' num2str(c)]);
        if strcmp(class(temp_rfr),'TreeBagger')
            fco2_rfr_tmp(:,c) = ...
                predict(temp_rfr,Vars_array.(region{n}).X_clust(:,Vars_array.(region{n}).idx_vars));
        else
            fco2_rfr_tmp(:,c) = NaN;
        end

        % reshape probabilities
        fco2_gmm_probs(:,c) = ...
            Clusts_grid.(region{n}).probabilities.(['c' num2str(c)])(Preds_grid.(region{n}).idx_clust);
    
    end
        
    % Clean up
    clear Clusts_grid Mods temp_rfr temp_nn c

    % average RFR result across clusters
    fco2_rfr = sum(fco2_rfr_tmp.*fco2_gmm_probs,2,'omitnan');

    % clean up
    clear fco2_rfr_tmp fco2_gmm_probs fco2_nn_temp fco2_nn

    %% re-grid
    % grid information
    OAI_grid.(region{n}).lon = Preds_grid.(region{n}).lon;
    OAI_grid.(region{n}).lat = Preds_grid.(region{n}).lat;
    OAI_grid.(region{n}).lim = Preds_grid.(region{n}).lim;
    OAI_grid.(region{n}).dim = Preds_grid.(region{n}).dim;
    OAI_grid.(region{n}).month = Preds_grid.(region{n}).month;
    OAI_grid.(region{n}).year = Preds_grid.(region{n}).year;
    OAI_grid.(region{n}).month_of_year = Preds_grid.(region{n}).month_of_year;
    OAI_grid.(region{n}).idxspc = Preds_grid.(region{n}).idxspc;

    % assemble fCO2 estimates on grid
    OAI_grid.(region{n}).fCO2 = nan(Preds_grid.(region{n}).dim.x,...
                        Preds_grid.(region{n}).dim.y,...
                        Preds_grid.(region{n}).dim.z);
    OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).idx_clust) = fco2_rfr;
    % remove fCO2 where ice covers more than 50% of grid cell
    OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).IceC > 0.5) = NaN;

    %% plot estimated fCO2
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).fCO2,...
        parula(20),'fCO2','Surface {\itf}CO_{2} (\muatm)',region{n});
    plot_regional_gif(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).fCO2,parula(20),'fCO2',...
        'Surface {\itf}CO_{2} (\muatm)',OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));

    %% save estimated fCO2 grid
    save(['Data/' region{n} '/ML_fCO2'],'OAI_grid','-v7.3');

    %% clean up
    clear archs fco2_rfr fco2_rfr_tmp fCO23Didx OAI_grid Preds_grid Vars_array

end

%% plot fCO2 across full region
plot_temporal_mean_full(295,475,10,parula(18),'fCO2','Sea Surface {\itf}CO_{2}',region,lme_shape,lme_idx)

%% plot fCO2 across full region (seasonally)
plot_temporal_mean_full_seas(295,475,10,parula(18),'fCO2','Sea Surface {\itf}CO_{2}',region,lme_shape,lme_idx)

%% plot gif of fCO2 across full region
plot_full_gif(295,475,parula(18),'fCO2','Sea Surface {\itf}CO_{2}',region,lme_shape,lme_idx);
