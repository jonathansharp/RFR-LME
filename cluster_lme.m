% This function uses the specified satellite and reanalysis variables to
% define time-varying clusters within which to train machine learning
% algorithms via a Gaussian mixture model method.

function cluster_lme(vrs,num_groups,region,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

for n = 1:length(region)

    % load variables for clustering
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    load(['Data/LME_Data/' vrs '_clustering_arrays_' region{n}],'LME_clustering');

    % normalize X to have mean of zero and st. dev. of 1
    X_norm = normalize(LME_clustering.monthly);
    
%     % determine covariance among cluster variables
%     idx = sum(~isnan(LME_clustering.(region{n})),2) == 2;
%     [R,P] = corrcoef(LME_clustering.(region{n})(idx,:));
%     R = [R(1,2),R(1,3),R(2,3)];
%     P = [P(1,2),P(1,3),P(2,3)];
%     R_idx = P < 0.05;
%     if sum(R_idx) >= 2
%         Sigma = 'full';
%     else
%         Sigma = 'diagonal';
%     end

    % fit GMM from predictor variables
    options = statset('MaxIter',1000,'Display','off');
    idx = ~any(isnan(X_norm),2);
    gmm = fitgmdist(X_norm(idx,:),num_groups(n),...
        'CovarianceType','full','SharedCovariance',false,...
        'Replicates',10,'Options',options,'RegularizationValue',0);
    [LME_GMM.clusters,~,LME_GMM.probs] = cluster(gmm,X_norm);
    
    % save BIC and AIC
    LME_GMM.BIC = gmm.BIC;
    LME_GMM.AIC = gmm.AIC;

%     % calculate silhouette score
%     sil = silhouette(X_norm,clusters);
    
%     % plot silhouette score
%     hold on
%     yL = ylim;
%     plot([mean(sil,'omitnan') mean(sil,'omitnan')],[yL(1) yL(2)],'k--','linewidth',2);
%     title([reg ', Clusters = ' num2str(num_groups)]);
%     export_fig(gcf,['Figures/gmm_validate_' reg '_' num2str(num_groups) '.png'],'-transparent')
%     close
    
    % fill group grid
    idx_clust = ~isnan(LME.SSS) & LME.idxspc;
    LME_GMM.group3D = nan(LME.dim.x,LME.dim.y,LME.dim.z);
    LME_GMM.group3D(idx_clust) = LME_GMM.clusters;

    % fill probability grids
    for g = 1:num_groups(n)
        LME_GMM.prob3D.(['c' num2str(g)]) = nan(LME.dim.x,LME.dim.y,LME.dim.z);
        LME_GMM.prob3D.(['c' num2str(g)])(idx_clust) = LME_GMM.probs(:,g);
    end

    % create animation
    if plot_option == 1
        create_animation([region{n} '_groups'],['GMM_' num2str(num_groups(n))],...
            datenum(LME.year,LME.month,15),LME.lat,LME.lon,...
            LME_GMM.group3D,jet,[1 num_groups(n)],'GMM Groups','');
    end

%     % fit GMM from predictor variable variability
%     options = statset('MaxIter',1000,'Display','off');
%     idx = ~any(isnan(X_norm),2);
%     gmm = fitgmdist(X_norm(idx,:),8,...
%         'CovarianceType','full',...
%         'SharedCovariance',false,'Replicates',1,...
%         'Options',options,'RegularizationValue',0);
%     [LME_GMM.clusters,~,LME_GMM.probs] = cluster(gmm,X_norm);

    % save variables for clustering and training in each LME
    save(['Data/LME_Data/' vrs '_GMM_' region{n}],'LME_GMM');
%    save(['Data/LME_Data/' vrs '_GMM_var_' region{n}],'LME_GMM_var');

end

end
