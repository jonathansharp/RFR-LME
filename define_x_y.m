% This function defines and assembles training (X) and target (pCO2)
% variables for  US large Marine Ecosystems in preparation for clustering
% via self-organizing maps and algorithm training.

function define_x_y(vrs,clust_vars,pred_vars,clust_vars_arc,...
    pred_vars_arc,clust_dims,pred_dims,region)

for n = 1:length(region)

    % load gridded pco2 and predictor data in individual LMEs
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');

    %% Compute regional average total and seasonal data coverage
    disp(['Processing data for clustering and algorithm training (' region{n} ')']);
    disp(['Total Coverage: ' sprintf('%.1f',mean(100.*LME.num_months(LME.idxspc(:,:,1))./LME.dim.z)) '%']);
    disp(['Seasonal Coverage: ' sprintf('%.1f',mean(100.*LME.num_months_clim(LME.idxspc(:,:,1))./12)) '%']);

    %% Create X and Y variables, as well as indices

    % create temporary versions of dimensions
    all_dims = unique([clust_dims pred_dims]);
    for v = 1:length([clust_dims pred_dims])
        if strcmp(all_dims{v},'lon')
            LME.([all_dims{v} '_grid']) = ...
                repmat(LME.(all_dims{v}),1,...
                    LME.dim.y,LME.dim.z);
        elseif strcmp(all_dims{v},'lat')
            LME.([all_dims{v} '_grid']) = ...
                repmat(LME.(all_dims{v})',...
                    LME.dim.x,1,LME.dim.z);
        elseif strcmp(all_dims{v},'sin_month_of_year')
            LME.([all_dims{v} '_grid']) = ...
                repmat(permute(sin((2.*pi.*LME.month_of_year)./12),[3 2 1]),...
                    LME.dim.x,LME.dim.y,1);
        elseif strcmp(all_dims{v},'cos_month_of_year')
            LME.([all_dims{v} '_grid']) = ...
                repmat(permute(cos((2.*pi.*LME.month_of_year)./12),[3 2 1]),...
                    LME.dim.x,LME.dim.y,1);
        elseif strcmp(all_dims{v},'dist') || strcmp(all_dims{v},'bathy')
            % skip
        else
            LME.([all_dims{v} '_grid']) = ...
                repmat(permute(LME.(all_dims{v}),[3 2 1]),...
                    LME.dim.x,LME.dim.y,1);
        end
    end

    % replicate bathymetry
    LME.Bathy = repmat(LME.Bathy,1,1,LME.dim.z);

    % calculate distance from shore
    LME.dist_grid = dist2coast(LME.lat_grid,LME.lon_grid);

    % create indices
    LME_clustering.idx = ~isnan(LME.SSS) & LME.idxspc; % within LME with valid SSS 
    LME_clustering.idx_var = ~isnan(LME.SSS(:,:,1)) & LME.idxspc(:,:,1); % within LME with valid SSS (2D)
    LME_training.idx = ~isnan(LME.fco2_ave_wtd) & LME.idxspc; % within LME with valid fCO2 observation
    LME_prediction.idx = ~isnan(LME.SSS) & LME.idxspc; % within LME with valid SSS 

    % calculate variability
    if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS')
        all_vars = unique([clust_vars_arc pred_vars_arc]);
    else
        all_vars = unique([clust_vars pred_vars]);
    end
    for v = 1:length(all_vars)
        LME.([pred_vars{v} '_var']) = ...
            std(LME.(pred_vars{v}),[],3);
    end

    % create X and Y for clustering
    if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS')
        clust_all = [clust_dims clust_vars_arc];
    else
        clust_all = [clust_dims clust_vars];
    end
    LME_clustering.headers = clust_all;
    LME_clustering.monthly = nan(sum(LME_clustering.idx(:)),length(clust_all));
    LME_clustering.variability = nan(sum(LME_clustering.idx_var(:)),length(clust_all));
    for v = 1:length(clust_all)
        if v <= length(clust_dims)
            LME_clustering.monthly(:,v) = LME.([clust_all{v} '_grid'])(LME_clustering.idx);
        else
            LME_clustering.monthly(:,v) = LME.(clust_all{v})(LME_clustering.idx);
        end
        if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS') 
            LME_clustering.variability(:,v) = LME.([clust_vars_arc{v} '_var'])(LME_clustering.idx_var);
        else
            LME_clustering.variability(:,v) = LME.([clust_vars{v} '_var'])(LME_clustering.idx_var);
        end
    end

    % create X and Y for algorithm training
    if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS') 
        pred_all = [pred_dims pred_vars_arc];
    else
        pred_all = [pred_dims pred_vars];
    end
    LME_training.headers = pred_all;
    LME_training.x = nan(sum(LME_training.idx(:)),length(pred_all));
    LME_training.y = nan(sum(LME_training.idx(:)),1);
    for v = 1:length(unique(pred_all))
        if v <= length(pred_dims)
            LME_training.x(:,v) = LME.([pred_all{v} '_grid'])(LME_training.idx);
        else
            LME_training.x(:,v) = LME.(pred_all{v})(LME_training.idx);
        end
    end
    LME_training.y = LME.fco2_ave_wtd(LME_training.idx);

    % create X and Y for prediction
    if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS')
        pred_all = [pred_dims pred_vars_arc];
    else
        pred_all = [pred_dims pred_vars];
    end
    LME_prediction.headers = pred_all;
    LME_prediction.x = nan(sum(LME_prediction.idx(:)),length(pred_all));
    for v = 1:length(unique(pred_all))
        if v <= length(pred_dims)
            LME_prediction.x(:,v) = LME.([pred_all{v} '_grid'])(LME_prediction.idx);
        else
            LME_prediction.x(:,v) = LME.(pred_all{v})(LME_prediction.idx);
        end
    end

    % save variables for clustering and training in each LME
    save(['Data/LME_Data/' vrs '_clustering_arrays_' region{n}],'LME_clustering');
    save(['Data/LME_Data/' vrs '_training_arrays_' region{n}],'LME_training');
    save(['Data/LME_Data/' vrs '_prediction_arrays_' region{n}],'LME_prediction');

end


end
