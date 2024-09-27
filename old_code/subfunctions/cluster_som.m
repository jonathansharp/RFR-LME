function group3D = cluster_som(X,idx_vars,lat,lon,time,idx_clust,dim1,dim2,type)

% Cluster data with self-organizing map

% Create a network
net = selforgmap([dim1 dim2]);

% Train the network
[net,~] = train(net,X(:,idx_vars)');
% Plot hits
figure; plotsomhits(net,X(:,idx_vars)');
% Test the network with input data
y = net(X(:,idx_vars)');

% Cluster input data into groups
group = nan(size(y,2),1);
for n=1:size(y,2)
    group(n) = find(y(:,n)==1);
end

% Plot histogram
figure; histogram(group);

% Preallocate group grid
if strcmp(type,'abs')
    group3D = nan(length(lon),length(lat),length(time));
elseif strcmp(type,'var')
    group3D = nan(length(lon),length(lat));
end

% Reshape groups to original data size
group3D = group3D(:);
group3Didx = find(idx_clust);
group3D(group3Didx) = group;
if strcmp(type,'abs')
    group3D = reshape(group3D,length(lon),length(lat),length(time));
elseif strcmp(type,'var')
    group3D = reshape(group3D,length(lon),length(lat));
end

clear Index SOMinputs net tr y group3Didx
