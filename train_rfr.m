% This function 

function train_rfr(vrs,num_groups,pred_dims,pred_vars,pred_vars_arc,...
    region,thresh,nTrees,minLeafSize,numpredictors)

for n = 1:length(region)

    % load variables for training
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    load(['Data/LME_Data/' vrs '_GMM_' region{n}],'LME_GMM');
    load(['Data/LME_Data/' vrs '_training_arrays_' region{n}],'LME_training');

    % define indices
    idx_mod = ~isnan(LME.fco2_ave_wtd) & LME.idxspc;
    LME_training_idx = false(sum(idx_mod(:)),num_groups(n));
    LME_probabilities = nan(sum(idx_mod(:)),num_groups(n));
    for c = 1:num_groups(n)
        LME_training_idx(:,c) = LME_GMM.prob3D.(['c' num2str(c)])(idx_mod) > thresh;
        LME_probabilities(:,c) = LME_GMM.prob3D.(['c' num2str(c)])(idx_mod);
    end

    % pre-allocate error stats
    rfr.avg_err = nan(num_groups(n),1);
    rfr.med_err = nan(num_groups(n),1);
    rfr.rmse = nan(num_groups(n),1);
    rfr.mae = nan(num_groups(n),1);
    rfr.r2 = nan(num_groups(n),1);
    
    % pre-allocate test pco2s and deltas
    rfr.Y_fit.all = nan(length(LME_training.y),num_groups(n));
    rfr.delta.all = nan(length(LME_training.y),num_groups(n));
    rfr.pred_imp.all = nan(length(LME_training.y),num_groups(n));

    % fit rfr
    for c = 1:num_groups(n) % for each cluster

    % cluster label
    clab = ['c' num2str(c)];

        if sum(LME_training_idx(:,c)) > 40 % need at least 40 observations to fit model (arbitrary)
            
            % define X and Y datasets for each cluster
            X_temp = LME_training.x(LME_training_idx(:,c),:);
            Y_temp = LME_training.y(LME_training_idx(:,c),:);
        
            % set up cross-validation
            numFolds = 5;
            sz = length(Y_temp);
            ints = randperm(length(Y_temp))';
        
            % pre-allocate
            rfr.Y_fit.(clab) = nan(size(Y_temp));
            rfr.delta.(clab) = nan(size(Y_temp));
            rfr.pred_imp.(clab) = nan(size(Y_temp));
        
            %% for each fold
            for f = 1:numFolds

                % fold label
                flab = ['f' num2str(f)];
                % index for cross-validation training and test
                idx_test = ints > (f-1) * (sz/numFolds) & ints <= f * (sz/numFolds);
                idx_train = ~idx_test;
                
                % get train and test data for this fold
                X_train = X_temp(idx_train,:);
                Y_train = Y_temp(idx_train);
                X_test = X_temp(idx_test,:);
                Y_test = Y_temp(idx_test);
        
                % construct random forest
                rfr_cv.(clab).(flab) = ...
                    TreeBagger(nTrees,X_train,Y_train,'Method','regression',...
                        'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors,...
                        'OOBPrediction','off','OOBPredictorImportance','off');
            
                % test random forest
                rfr.Y_fit.(clab)(idx_test) = ...
                    predict(rfr_cv.(clab).(flab),X_test);
                rfr.delta.(clab)(idx_test) = rfr.Y_fit.(clab)(idx_test)-Y_test;
            
            end
        
            % save network error statistics for each fold
            rfr.avg_err(c) = mean(rfr.delta.(clab));
            rfr.med_err(c) = median(rfr.delta.(clab));
            rfr.rmse(c) = sqrt(mean(rfr.delta.(clab).^2));
            rfr.mae(c) = median(abs(rfr.delta.(clab)));
            rfr.r2(c) = corr(rfr.Y_fit.(clab),Y_temp).^2;
        
            % save fitted values and errors across clusters
            rfr.Y_fit.all(LME_training_idx(:,c),c) = rfr.Y_fit.(clab);
            rfr.delta.all(LME_training_idx(:,c),c) = rfr.delta.(clab);
        
            % construct random forest with all data
            rfr.(clab) = ...
                TreeBagger(nTrees,X_temp,Y_temp,'Method','regression',...
                    'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors,...
                    'OOBPrediction','on','OOBPredictorImportance','on');
        
            % Plot Out-of-Bag MSE based on tree number
            figure('visible','on');
            plot(sqrt(oobError(rfr.(clab))),'k','LineWidth',2);
            xlabel('Number of Grown Trees');
            ylabel('Out-of-Bag Root Mean Squared Error');
            % save figure
            if ~isfolder(['Figures/rfr_validate/' region{n}])
                mkdir(['Figures/rfr_validate/' region{n}]); end
            export_fig(gcf,['Figures/rfr_validate/' region{n} '/' clab '_OOB_RMSE.png'],'-transparent');
            close
            
            % Plot importance of each predictor
            figure('visible','on');
            set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
            set(gca,'fontsize',22);
            bar(rfr.(clab).OOBPermutedPredictorDeltaError,'k');
            set(gca,'fontsize',16);
            xlabel('Predictors');
            ylabel('Out-of-Bag Feature Importance');
            if strcmp(region{n},'BS') || strcmp(region{n},'NBCS') || strcmp(region{n},'EBS') 
                xticks(1:length([pred_dims pred_vars_arc]));
                xticklabels([pred_dims pred_vars_arc]);
            else
                xticks(1:length([pred_dims pred_vars]));
                xticklabels([pred_dims pred_vars]);
            end
            % export figure
            export_fig(gcf,['Figures/rfr_validate/' region{n} '/' clab '_Pred_Imp.png'],'-transparent');
            close
    
        else
    
            rfr.Y_fit.(clab) = nan;
            rfr.delta.(clab) = nan;
            rfr.pred_imp.(clab) = nan;
            rfr.avg_err(c) = nan;
            rfr.med_err(c) = nan;
            rfr.rmse(c) = nan;
            rfr.mae(c) = nan;
            rfr.r2(c) = nan;
            rfr.(clab) = nan;
    
        end

    end

% average error statistics across clusters
LME_probabilities(~LME_training_idx) = NaN;
rfr.Y_fit.all(:,c+1) = sum(rfr.Y_fit.all.*LME_probabilities,2,'omitnan')./...
    sum(LME_probabilities,2,'omitnan');
rfr.delta.all(:,c+1) = sum(rfr.delta.all.*LME_probabilities,2,'omitnan')./...
    sum(LME_probabilities,2,'omitnan');

% rfr error statistics across all clusters
idx = ~isnan(rfr.Y_fit.all(:,c+1));
rfr.avg_err(c+1) = mean(rfr.delta.all(:,c+1),'omitnan');
rfr.med_err(c+1) = median(rfr.delta.all(:,c+1),'omitnan');
rfr.rmse(c+1) = sqrt(mean(rfr.delta.all(:,c+1).^2,'omitnan'));
rfr.mae(c+1) = median(abs(rfr.delta.all(:,c+1)),'omitnan');
rfr.r2(c+1) = corr(rfr.Y_fit.all(idx,c+1),LME_training.y(idx)).^2;

% place errors on grid
rfr.delta_grid = nan(size(LME_training.idx));
rfr.delta_grid_abs = nan(size(LME_training.idx));
rfr.delta_grid(LME_training.idx) = rfr.delta.all(:,end);
rfr.delta_grid_abs(LME_training.idx) = abs(rfr.delta.all(:,end));
% add dimensions to grid
rfr.lat = LME.lat;
rfr.lon = LME.lon;
rfr.month = LME.month;

%% plot 2D histogram of delta pco2 for RFRs
x_edges = 0:5:1000;
x_mids = 2.5:5:997.5;
y_edges = -1000:10:1000;
y_mids = -995:10:995;
plot_delta(x_edges,x_mids,y_edges,y_mids,rfr.Y_fit.all(:,end),...
    rfr.delta.all(:,end),region{n},[200 800],[-300 300],'RFR');

%% plot 2D histogram of pco2 vs. pco2 for RFRs
x_edges = 0:5:1000;
x_mids = 2.5:5:997.5;
y_edges = 0:5:1000;
y_mids = 0:5:1000;
plot_mod_v_meas(x_edges,x_mids,y_edges,y_mids,LME_training.y,...
    rfr.Y_fit.all(:,end), rfr.delta.all(:,end),...
    region{n},[200 600],[200 600],'RFR');

% save rfr model and error statistics
if ~isfolder(['Models/' region{n}]); mkdir(['Models/' region{n}]); end
save(['Models/' region{n} '/us_lme_models_' vrs '_c' num2str(c)],'rfr');

% display information
disp(['RFR trained for ' region{n}]);

end

end
