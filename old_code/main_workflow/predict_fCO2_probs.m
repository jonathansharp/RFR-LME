% Predict fCO2 on grid
% 
% This script uses neural networks and random forest regressions to predict
% fCO2 in time-varying clusters defined by a self-organizing map method for
% US large Marine Ecosystems.
% 
% Written by J.D. Sharp: 9/15/22
% Last updated by J.D. Sharp: 10/17/24
% 

% Set GMM parameters
set_gmm_options

% include moorings?
if include_moorings == 0; exts = '_no_moorings'; else; exts = ''; end

for n = 7%1:length(region)

    %% display status
    disp(['Predicting fCO2 (' region{n} ')']);

    %% load gridded fCO2 and predictors
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/variable_arrays' exts],'Vars_array');
    load(['Data/' region{n} '/gridded_clusters'],'Clusts_grid');
    load(['Data/' region{n} '/us_lme_models' exts],'Mods');

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

    %% process ice-covered grid cells
    % remove fCO2 where ice covers more than 50% of grid cell
    OAI_grid.(region{n}).idxspc(Preds_grid.(region{n}).IceC > 0.5) = 0;
    OAI_grid.(region{n}).fCO2(~OAI_grid.(region{n}).idxspc) = NaN;
%     % fill removed gaps
%     if any((Preds_grid.(region{n}).IceC(:) > 0.5))
%         for a = 1:OAI_grid.(region{n}).dim.x
%             for b = 1:OAI_grid.(region{n}).dim.y
%                 if any((Preds_grid.(region{n}).IceC(a,b,:) > 0.5))
%                     % define fCO2 and time
%                     y = squeeze(OAI_grid.(region{n}).fCO2(a,b,:));
%                     t = (1:OAI_grid.(region{n}).dim.z)';
%                     % index to where fCO2 is NaN
%                     idx = ~isnan(y);
%                     % find distances between NaNs and values
%                     dists = nan(size(y));
%                     for d = 1:length(dists)
%                         dists(d) = min(abs(d-find(idx)));
%                     end
%                     % interpolate using pchip
%                     %yt = interp1(t(idx),y(idx),t,'pchip');
%                     % interpolate using nearest
%                     %yt = interp1(t(idx),y(idx),t,'nearest');
%                     % interpolate using pchip first, then nearest
%                     idx_pchip = dists <= 2;
%                     y_pchip = interp1(t(idx),y(idx),t(idx_pchip),'pchip');
%                     yt = interp1(t(idx_pchip),y_pchip,t,'nearest','extrap');
%                     % plot
%                     %figure; set(gcf,'Position',[400 400 1120 420]); plot(t,nan(size(t)),t,y);
%                     %figure; set(gcf,'Position',[400 400 1120 420]); plot(t,yt,t,y);
%                     % add interpolated fCO2 to grid
%                     OAI_grid.(region{n}).fCO2(a,b,:) = yt;
%                 end
%             end
%         end
%     end
%     clear a b y t idx dists idx_pchip y_pchip idx_nearest yt

    %% convert fCO2 to pCO2
    % calculate fugacity factor
    TempK = Preds_grid.(region{n}).SST + 273.15;
    Delta = 57.7 - 0.118.*TempK;
    b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
    RGasConstant = 83.14462618;
    FugFac = exp((b + 2.*Delta).*Preds_grid.(region{n}).mslp./(RGasConstant.*TempK));
    OAI_grid.(region{n}).pCO2 = OAI_grid.(region{n}).fCO2./FugFac;
    clear TempK Delta b RGasConstant FugFac

    %% save estimated fCO2 grid
    save(['Data/' region{n} '/ML_fCO2' exts],'OAI_grid','-v7.3');

    %% clean up
    clear archs fco2_rfr fco2_rfr_tmp fCO23Didx OAI_grid Preds_grid Vars_array

end
