function group3D = cluster_data(X,idx_vars,lat,lon,time,idx_clust,reg)

% Cluster data with self-organizing map

% Create a network
dimension1 = 3;
dimension2 = 3;
net = selforgmap([dimension1 dimension2]);

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

% Reshape groups to original data size
group3D = nan(length(lon),length(lat),length(time));
group3D = group3D(:);
group3Didx = find(idx_clust);
group3D(group3Didx) = group;
group3D = reshape(group3D,length(lon),length(lat),length(time));

clear Index SOMinputs net tr y group3Didx

% Plot groups on map
figure;
worldmap([min(lat) max(lat)],[min(lon) max(lon)]);
pcolorm(lat,lon,mode(group3D,3)');
colormap(jet(9));
title('Most frequent group from self-organizing map method');
c=colorbar;
caxis([0.5 9.5]);
c.Label.String = 'Groups';
c.TickLength = 0;

% save figure
exportgraphics(gcf,['Figures/' reg '_SOM.png']);
close
