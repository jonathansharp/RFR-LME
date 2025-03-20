% This function uses the specified satellite and reanalysis variables to
% define time-varying clusters within which to train machine learning
% algorithms via a Gaussian mixture model method.

function cluster_lme(vrs,num_groups,region)

% load variables for clustering
load(['Data/' vrs '_clustering_array']);

for n = 1:length(region)

    for en = 1:size(num_groups,2)

%         %% load gridded pCO2 and predictors
%         load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
%         load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    
        %% cluster gridded data
        disp(['Clustering data for ' region{n}]);
        
%         % define variable index, but replace Chlorophyll with wind speed for
%         % Bering-Chukchi and Beaufort Seas
%         if strcmp(region{n},'BS') || strcmp(region{n},'NBCS')
%             vars = {'MSLP' 'SST' 'Wind'};
%             % vars = {'MSLP' 'SST' 'Wind' 'Lat.' 'Lon.'};
%         else
%             vars = {'MSLP' 'SST' 'CHL'};
%             % vars = {'MSLP' 'SST' 'CHL' 'Lat.' 'Lon.'};
%         end

%         % determine covariance among cluster variables
%         idx = sum(~isnan(LME_clustering.(region{n})),2) == 2;
%         [R,P] = corrcoef(LME_clustering.(region{n})(idx,:));
%         R = [R(1,2),R(1,3),R(2,3)];
%         P = [P(1,2),P(1,3),P(2,3)];
%         R_idx = P < 0.05;
%         if sum(R_idx) >= 2
%             Sigma = 'full';
%         else
%             Sigma = 'diagonal';
%         end

        % plotmatrix(clust_temp(:,idx_vars));
        % obtain gridded group numbers
        rng(3); % for reproducibility
        % Cluster observations
        done = 0; % set done to zero
        counter = 1; % set counter to 1
        while done == 0 && counter < 100
            try
                counter = counter+1; % add attempt to counter
                [Clusts_grid.(region{n}).groups,...
                    Clusts_grid.(region{n}).probabilities,...
                    Clusts_grid.(region{n}).BIC,...
                    Clusts_grid.(region{n}).SIL] = ...
                    cluster_gmm(clust_temp,idx_vars,...
                        Preds_grid.(region{n}).lat,Preds_grid.(region{n}).lon,...
                        Preds_grid.(region{n}).month_of_year,...
                        idx_clust_temp,...
                        num_groups(n,en),strat,options,Sigma,...
                        SharedCovariance,RegularizationValue,region{n});
                done = 1; % successful cluster
            catch
                rng(randi(100)) % reset random seed and try again if there's an error
            end
        end
        levs = num_groups(n,en);
        type = 'gmm';

        % replicate groups and probabilities over time
        if strcmp(strat,'var')
            Clusts_grid.(region{n}).groups = ...
                repmat(Clusts_grid.(region{n}).groups,1,1,...
                Preds_grid.(region{n}).dim.z);
            for c = 1:num_groups(n,en)
                Clusts_grid.(region{n}).probabilities.(['c' num2str(c)]) = ...
                    repmat(Clusts_grid.(region{n}).probabilities.(['c' num2str(c)]),1,1,...
                    Preds_grid.(region{n}).dim.z);
            end
        end
    
        % clean up
        clear vars idx_vars idx_clust clust_temp idx_clust_temp

        %% Save gridded cluster data for all LMEs
        if gmm_test_idx == 0
            % save(['Data/' region{n} '/gridded_clusters_' num2str(en)],'Clusts_grid','-v7.3');
            save(['Data/' region{n} '/gridded_clusters'],'Clusts_grid','-v7.3');
        elseif gmm_test_idx == 1
            save(['Data/' region{n} '/gridded_clusters_test'],'Clusts_grid','-v7.3');
        end
    
        % clean up
        clear levs type c Preds_grid Vars_array Clusts_grid

    end

end

end
