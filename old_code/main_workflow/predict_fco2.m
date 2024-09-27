% Predict fCO2 on grid
% 
% This script uses neural networks and random forest regressions to predict
% fCO2 in time-varying clusters defined by a self-organizing map method for
% US large Marine Ecosystems.
% 
% Written by J.D. Sharp: 9/15/22
% Last updated by J.D. Sharp: 1/23/22
% 

% this script defines the bounds of the eighteen LMEs
define_regions

for n = 1:length(region)

    %% display status
    disp(['Predicting fCO2 (' region{n} ')']);

    %% load gridded fCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');

    %% pre-allocate
    % archs = length(fields(Mods.(region{n}).nn.c1));
    % fco2_nn_temp = nan(length(Preds_grid.(region{n}).idx_group_clust),archs);
    fco2_rfr_tmp = nan(length(Vars_array.(region{n}).X_clust),size(num_groups,2));

    % for each set of clusters
    for en = 1:size(num_groups,2)

        %% load gridded models and clusters
        load(['Data/' region{n} '/gridded_clusters_' num2str(en)],'Clusts_grid');
        load(['Data/' region{n} '/us_lme_models_' num2str(en)],'Mods');
    
        %% define group indices
        Clusts_grid.(region{n}).idx_group_clust = ...
            Clusts_grid.(region{n}).groups(Preds_grid.(region{n}).idx_clust);
    
        %% predict fCO2 for each cluster
        for c = 1:max(Clusts_grid.(region{n}).idx_group_clust)
    
            % cluster index
            idx_tmp = Clusts_grid.(region{n}).idx_group_clust == c;
    
    %         % apply neural networks
    %         for a = 1:archs
    %             temp_nn = ...
    %                 Mods.(region{n}).nn.(['c' num2str(c)]).(['a' num2str(a)]);
    %             fco2_nn_temp(idx_tmp,a) = ...
    %                 temp_nn(Vars_array.(region{n}).X_clust(idx_tmp,Vars_array.(region{n}).idx_vars)')';
    %         end
    
            % apply random forest regressions
            temp_rfr = Mods.(region{n}).rfr.(['c' num2str(c)]);
            if sum(idx_tmp > 0)
                if strcmp(class(temp_rfr),'TreeBagger')
                    fco2_rfr_tmp(idx_tmp,en) = ...
                        predict(temp_rfr,Vars_array.(region{n}).X_clust(idx_tmp,...
                        Vars_array.(region{n}).idx_vars));
                else
                    fco2_rfr_tmp(idx_tmp,en) = NaN;
                end
            end
        
        end
        
        % Clean up
        clear Clusts_grid Mods

    end

    % average RFR result across clusters
    fco2_rfr = mean(fco2_rfr_tmp,2,'omitnan');

    % average neural network results across architectures
%     fco2_nn = mean(fco2_nn_temp,2);

    % average NN and RFR result
%     fco2_avg = mean([fco2_nn fco2_rfr],2);

    % clean up
    clear c idx_tmp temp_nn fco2_nn_temp fco2_nn temp_rfr

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
%     fCO2_grid.(region{n}).fCO2 = fCO2_grid.(region{n}).fCO2(:);
%     fCO23Didx = find(Preds_grid.(region{n}).idx_clust);
    OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).idx_clust) = fco2_rfr;
%     OAI_grid.(region{n}).fCO2 = ...
%         reshape(OAI_grid.(region{n}).fCO2,OAI_grid.(region{n}).dim.x,...
%                               OAI_grid.(region{n}).dim.y,...
%                               OAI_grid.(region{n}).dim.z);

%     %% smooth in two dimensions
%     load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
%     OAI_grid.(region{n}).fCO2_sm = nan(size(OAI_grid.(region{n}).fCO2));
%     for t = 1:OAI_grid.(region{n}).dim.z
%         OAI_grid.(region{n}).fCO2_sm(:,:,t) = ...
%             movmean(OAI_grid.(region{n}).fCO2(:,:,t),3,1);
%         OAI_grid.(region{n}).fCO2_sm(:,:,t) = ...
%             movmean(OAI_grid.(region{n}).fCO2_sm(:,:,t),3,2);
%     end

    % average fCO2 across

    %% plot estimated fCO2
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).fCO2,...
        parula(20),'fCO2','Surface {\itf}CO_{2} (\muatm)',region{n});

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
