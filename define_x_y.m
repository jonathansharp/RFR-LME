% This function defines and assembles training (X) and target (pCO2)
% variables for  US large Marine Ecosystems in preparation for clustering
% via self-organizing maps and algorithm training.

function define_x_y(vrs,clust_vars,pred_vars,clust_dims,pred_dims,region)

% load gridded pco2 and predictor data in individual LMEs
load(['Data/' vrs '_gridded_lme'],'LME_grid');

for n = 1:length(region)

    %% Compute regional average total and seasonal data coverage
    disp(['Processing data for clustering and algorithm training (' region{n} ')']);
    disp(['Total Coverage: ' sprintf('%.1f',mean(100.*LME_grid.(region{n}).num_months(LME_grid.(region{n}).idxspc(:,:,1))./LME_grid.(region{n}).dim.z)) '%']);
    disp(['Seasonal Coverage: ' sprintf('%.1f',mean(100.*LME_grid.(region{n}).num_months_clim(LME_grid.(region{n}).idxspc(:,:,1))./12)) '%']);

    %% Create X and Y variables, as well as indices

    % create temporary versions of dimensions
    all_dims = unique([clust_dims pred_dims]);
    for v = 1:length([clust_dims pred_dims])
        if strcmp(all_dims{v},'lon')
            LME_grid.(region{n}).([all_dims{v} '_grid']) = ...
                repmat(LME_grid.(region{n}).(all_dims{v}),1,...
                    LME_grid.(region{n}).dim.y,LME_grid.(region{n}).dim.z);
        elseif strcmp(all_dims{v},'lat')
            LME_grid.(region{n}).([all_dims{v} '_grid']) = ...
                repmat(LME_grid.(region{n}).(all_dims{v})',...
                    LME_grid.(region{n}).dim.x,1,LME_grid.(region{n}).dim.z);
        elseif strcmp(all_dims{v},'sin_month_of_year')
            LME_grid.(region{n}).([all_dims{v} '_grid']) = ...
                repmat(permute(sin((2.*pi.*LME_grid.(region{n}).month_of_year)./12),[3 2 1]),...
                    LME_grid.(region{n}).dim.x,LME_grid.(region{n}).dim.y,1);
        elseif strcmp(all_dims{v},'cos_month_of_year')
            LME_grid.(region{n}).([all_dims{v} '_grid']) = ...
                repmat(permute(cos((2.*pi.*LME_grid.(region{n}).month_of_year)./12),[3 2 1]),...
                    LME_grid.(region{n}).dim.x,LME_grid.(region{n}).dim.y,1);
        elseif strcmp(all_dims{v},'dist') || strcmp(all_dims{v},'bathy')
            % skip
        else
            LME_grid.(region{n}).([all_dims{v} '_grid']) = ...
                repmat(permute(LME_grid.(region{n}).(all_dims{v}),[3 2 1]),...
                    LME_grid.(region{n}).dim.x,LME_grid.(region{n}).dim.y,1);
        end
    end

    % calculate distance from shore
    LME_grid.(region{n}).dist = dist2coast(LME_grid.(region{n}).lat_grid,LME_grid.(region{n}).lon_grid);
    % calculate bathymetry
    LME_grid.(region{n}).bathy = bottom_depth(LME_grid.(region{n}).lat_grid,LME_grid.(region{n}).lon_grid);

    % create indices
    idx_clust = ~isnan(LME_grid.(region{n}).SSS) & ...
           repmat(LME_grid.(region{n}).idxspc(:,:,1),1,1,LME_grid.(region{n}).dim.z);
    idx_clust_var = ~isnan(LME_grid.(region{n}).SSS(:,:,1)) & ...
           LME_grid.(region{n}).idxspc(:,:,1);
    idx_mod = ~isnan(LME_grid.(region{n}).fco2_ave_wtd) & ...
           repmat(LME_grid.(region{n}).idxspc(:,:,1),1,1,LME_grid.(region{n}).dim.z);

    % calculate variability
    all_vars = unique([clust_vars pred_vars]);
    for v = 1:length(all_vars)
        LME_grid.(region{n}).([pred_vars{v} '_var']) = ...
            std(LME_grid.(region{n}).(pred_vars{v}),[],3);
    end

    % create X and Y for clustering
    LME_clustering.(region{n}) = nan(sum(idx_clust(:)),length(clust_vars));
    LME_clustering_var.(region{n}) = nan(sum(idx_clust_var(:)),length(clust_vars));
    for v = 1:length(clust_vars)
        LME_clustering.(region{n})(:,v) = LME_grid.(region{n}).(clust_vars{v})(idx_clust);
        LME_clustering_var.(region{n})(:,v) = LME_grid.(region{n}).([clust_vars{v} '_var'])(idx_clust_var);
    end

    % create X and Y for algorithm training
    LME_training.(region{n}) = nan(sum(idx_mod(:)),length(pred_vars));
    for v = 1:length(unique(pred_vars))
        LME_training.(region{n})(:,v) = LME_grid.(region{n}).(pred_vars{v})(idx_mod);
    end

    % add indices to gridded predictors
    LME_grid.(region{n}).idx_clust = idx_clust;
    LME_grid.(region{n}).idx_clust_var = idx_clust_var;
    LME_grid.(region{n}).idx_mod = idx_mod;

end

% save gridded pco2 and predictor data in individual LMEs
save(['Data/' vrs '_gridded_lme_processed'],'LME_grid','-v7.3');
save(['Data/' vrs '_clustering_array'],'LME_clustering');
save(['Data/' vrs '_clustering_array_var'],'LME_clustering_var');
save(['Data/' vrs '_training_array'],'LME_training');

end
