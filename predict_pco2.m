% Predict pCO2 on grid
% 
% This script uses neural networks and random forest regressions to predict
% pCO2 in time-varying clusters defined by a self-organizing map method for
% five US large Marine Ecosystems (Alaska, California Current, Insular
% Pacific / Hawaii, Gulf of Mexico / Caribbean, and US East Coast).
% 
% Written by J.D. Sharp: 9/15/22
% Last updated by J.D. Sharp: 9/28/22
% 

% define regions of interest
region = {'CCS' 'AK' 'EastCoast' 'GoM_Car' 'Hawaii'};

for n = 1:length(region)

    %% load gridded pCO2, predictors, and models
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    load(['Data/' region{n} '/us_lme_models'],'Mods');

    %% define group indices
    Preds_grid.(region{n}).idx_group_clust = ...
        Preds_grid.(region{n}).groups(Preds_grid.(region{n}).idx_clust);

    %% pre-allocate
    archs = length(fields(Mods.(region{n}).nn.c1));
    pco2_nn_temp = nan(length(Preds_grid.(region{n}).idx_group_clust),archs);
    pco2_rfr = nan(size(Preds_grid.(region{n}).idx_group_clust));

    %% predict pCO2 for each cluster
    disp(['Predicting pCO2 (' region{n} ')']);
    for c = 1:max(Preds_grid.(region{n}).idx_group_clust)

        % cluster index
        idx_tmp = Preds_grid.(region{n}).idx_group_clust == c;

        % apply neural networks
        for a = 1:archs
            temp_nn = ...
                Mods.(region{n}).nn.(['c' num2str(c)]).(['a' num2str(a)]);
            pco2_nn_temp(idx_tmp,a) = ...
                temp_nn(Vars_array.(region{n}).X_clust(idx_tmp,Vars_array.(region{n}).idx_vars)')';
        end

        % apply random forest regressions
        temp_rfr = Mods.(region{n}).rfr.(['c' num2str(c)]);
        pco2_rfr(idx_tmp) = predict(temp_rfr,Vars_array.(region{n}).X_clust(idx_tmp,Vars_array.(region{n}).idx_vars));
    
    end

    % average neural network results across architectures
    pco2_nn = mean(pco2_nn_temp,2);

    % average NN and RFR result
    pco2_avg = mean([pco2_nn pco2_rfr],2);

    %% re-grid
    % grid information
    OAI_grid.(region{n}).lon = Preds_grid.(region{n}).lon;
    OAI_grid.(region{n}).lat = Preds_grid.(region{n}).lat;
    OAI_grid.(region{n}).lim = Preds_grid.(region{n}).lim;
    OAI_grid.(region{n}).dim = Preds_grid.(region{n}).dim;
    OAI_grid.(region{n}).month = Preds_grid.(region{n}).month;
    OAI_grid.(region{n}).year = Preds_grid.(region{n}).year;
    OAI_grid.(region{n}).month_of_year = Preds_grid.(region{n}).month_of_year;
    % assemble pCO2 estimates on grid
    OAI_grid.(region{n}).pCO2 = nan(OAI_grid.(region{n}).dim.x,...
                        OAI_grid.(region{n}).dim.y,...
                        OAI_grid.(region{n}).dim.z);
    OAI_grid.(region{n}).pCO2 = OAI_grid.(region{n}).pCO2(:);
    pCO23Didx = find(Preds_grid.(region{n}).idx_clust);
    OAI_grid.(region{n}).pCO2(pCO23Didx) = pco2_avg;
    OAI_grid.(region{n}).pCO2 = ...
        reshape(OAI_grid.(region{n}).pCO2,OAI_grid.(region{n}).dim.x,...
                              OAI_grid.(region{n}).dim.y,...
                              OAI_grid.(region{n}).dim.z);

    % plot estimated pCO2
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).pCO2,...
        200,500,31,'pCO2_surf','Surface pCO_{2} (\muatm)',region{n});


    % save estimated fCO2 grid
    save(['Data/' region{n} '/ML_fCO2'],'OAI_grid','-v7.3');
    clearvars -except region

end
