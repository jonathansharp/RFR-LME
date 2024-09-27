function [group3D,sil] = cluster_kmeans(X,idx_vars,lat,lon,time,idx_clust,...
    num_groups,type,options)

% Cluster data with k-means clustering
[clusters,~] = kmeans(X(:,idx_vars),num_groups,'Options',options);

% Preallocate group grid
if strcmp(type,'abs')
    group3D = nan(length(lon),length(lat),length(time));
    sil = nan;
elseif strcmp(type,'var')
    group3D = nan(length(lon),length(lat));
    % calculate silhouette score
    s = silhouette(X(:,idx_vars),clusters);
    sil = mean(s);
end

% Reshape groups to original data size
group3D = group3D(:);
group3Didx = find(idx_clust);
group3D(group3Didx) = clusters;
if strcmp(type,'abs')
    group3D = reshape(group3D,length(lon),length(lat),length(time));
elseif strcmp(type,'var')
    group3D = reshape(group3D,length(lon),length(lat));
end

clear Index SOMinputs net tr y group3Didx
