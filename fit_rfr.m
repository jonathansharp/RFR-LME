% Fit Random Forest Regression
% 
% Written by J.D. Sharp: 9/2/22
% Last updated by J.D. Sharp: 9/15/22
% 

function [rfr,err,rmse,mae,r2,Y_fit,delta,pred_imp] = fit_rfr(X_mod,Y_mod,idx_vars,...
    idx_group,headers,nTrees,minLeafSize,numpredictors,region)

% pre-allocate error stats
err = nan(max(idx_group),1);
rmse = nan(max(idx_group),1);
mae = nan(max(idx_group),1);
r2 = nan(max(idx_group),1);

% pre-allocate test pco2s and deltas
Y_fit.all = nan(size(Y_mod));
delta.all = nan(size(Y_mod));
pred_imp.all = nan(size(Y_mod));

% for each cluster
for c = 1:max(idx_group)

    %% prepare

    % cluster label
    clab = ['c' num2str(c)];

    % cluster index
    idx_clust = idx_group==c;

    if sum(idx_clust) > 40 % need at least 40 observations to fit model (arbitrary)
        
        % define X and Y datasets for each cluster
        X_temp = X_mod(idx_clust,idx_vars);
        Y_temp = Y_mod(idx_clust);
        headers_temp = headers(idx_vars);
    
        % set up cross-validation
        numFolds = 5;
        sz = length(Y_temp);
        ints = randperm(length(Y_temp))';
    
        % pre-allocate
        Y_fit.(clab) = nan(size(Y_temp));
        delta.(clab) = nan(size(Y_temp));
        pred_imp.(clab) = nan(size(Y_temp));
    
    %     %% determine predictors to use
    %     idx_vars = 1:length(idx_vars);
    %     while length(idx_vars) > 10
    %         % get train data
    %         X_train = X_temp(:,idx_vars);
    %         % construct random forest with all data
    %         rfr_preds = ...
    %             TreeBagger(nTrees,X_train,Y_temp,'Method','regression',...
    %                 'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors,...
    %                 'OOBPrediction','on','OOBPredictorImportance','on');
    %         % determine least important predictor
    %         idx_vars(rfr_preds.OOBPermutedPredictorDeltaError == ...
    %             min(rfr_preds.OOBPermutedPredictorDeltaError)) = [];
    %     end
    
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
            try
            rfr_cv.(clab).(flab) = ...
                TreeBagger(nTrees,X_train,Y_train,'Method','regression',...
                    'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors,...
                    'OOBPrediction','off','OOBPredictorImportance','off');
            catch
                keyboard
            end
        
            % test random forest
            Y_fit.(clab)(idx_test) = ...
                predict(rfr_cv.(clab).(flab),X_test);
            delta.(clab)(idx_test) = Y_fit.(clab)(idx_test)-Y_test;
        
        end
    
        %% save network error statistics for each fold
        err(c) = mean(delta.(clab));
        rmse(c) = sqrt(mean(delta.(clab).^2));
        mae(c) = median(abs(delta.(clab)));
        r2(c) = corr(Y_fit.(clab),Y_temp).^2;
    
        %% save fitted values and errors across clusters
        Y_fit.all(idx_clust) = Y_fit.(clab);
        delta.all(idx_clust) = delta.(clab);
    
        %% construct random forest with all data
        rfr.(clab) = ...
            TreeBagger(nTrees,X_temp,Y_temp,'Method','regression',...
                'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors,...
                'OOBPrediction','on','OOBPredictorImportance','on');
    
        %% Plot Out-of-Bag MSE based on tree number
        figure('visible','off');
        plot(sqrt(oobError(rfr.(clab))),'k','LineWidth',2);
        xlabel('Number of Grown Trees');
        ylabel('Out-of-Bag Root Mean Squared Error');
        % save figure
        if ~isfolder('Figures/rfr_validate'); mkdir('Figures/rfr_validate'); end
        exportgraphics(gcf,['Figures/rfr_validate/' region '_' clab '_OOB_RMSE.png']);
        close
        
        %% Plot importance of each predictor
        figure('visible','off');
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gca,'fontsize',22);
        bar(rfr.(clab).OOBPermutedPredictorDeltaError,'k');
        set(gca,'fontsize',16);
        xlabel('Predictors');
        ylabel('Out-of-Bag Feature Importance');
        xticks(1:length(headers_temp));
        xticklabels(headers_temp);
        % export figure
        if ~isfolder('Figures/rfr_validate'); mkdir('Figures/rfr_validate'); end
        exportgraphics(gcf,['Figures/rfr_validate/' region '_' clab '_Pred_Imp.png']);
        close

    else

        Y_fit.(clab) = nan;
        delta.(clab) = nan;
        pred_imp.(clab) = nan;
        err(c) = nan;
        rmse(c) = nan;
        mae(c) = nan;
        r2(c) = nan;
        rfr.(clab) = nan;

    end

end

% save network error statistics across all clusters
err(c+1) = mean(delta.all);
rmse(c+1) = sqrt(mean(delta.all.^2));
mae(c+1) = median(abs(delta.all));
r2(c+1) = corr(Y_fit.all,Y_mod).^2;

% end function
end
