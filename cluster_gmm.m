function [group3D,prob3D,BIC] = cluster_gmm(X,idx_vars,lat,lon,time,idx_clust,...
    num_groups,type,options,Sigma,SharedCovariance,RegularizationValue)

% normalize X to have mean of zero and st. dev. of 1
X_norm = normalize(X(:,idx_vars));

% Fit GMM from predictor variables
gmm = fitgmdist(X_norm,num_groups,...
    'CovarianceType',Sigma,...
    'SharedCovariance',SharedCovariance,'Replicates',20,...
    'Options',options,'RegularizationValue',RegularizationValue);
[clusters,~,p] = cluster(gmm,X_norm);
BIC = gmm.BIC;
% calculate silhouette score
s = silhouette(X(:,idx_vars),clusters);
sil = mean(s);
% plot silhouette score


% Fill group grids
if strcmp(type,'abs')
    group3D = nan(length(lon),length(lat),length(time));
    group3D(idx_clust) = clusters;
elseif strcmp(type,'var')
    group3D = nan(length(lon),length(lat));
    group3D(idx_clust) = clusters;
end

% Fill probability grids
for g = 1:num_groups
    if strcmp(type,'abs')
        prob3D.(['c' num2str(g)]) = nan(length(lon),length(lat),length(time));
        prob3D.(['c' num2str(g)])(idx_clust) = p(:,g);
    elseif strcmp(type,'var')
        prob3D.(['c' num2str(g)]) = nan(length(lon),length(lat));
        prob3D.(['c' num2str(g)])(idx_clust) = p(:,g);
    end
end

clear Index SOMinputs net tr y group3Didx