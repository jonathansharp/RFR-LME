%% Cluster Analysis of means
% load and process data
data=readtable('IndsAndStats/MeansTable-15-Aug-2023.xls');
LMEs = data{:,1};
Vars = data.Properties.VariableNames;
data=table2array(data(:,[2 3 5 6 7 11]));
% Standardize data
nr = size(data,1);
data_stnd = (data - repmat(mean(data),nr,1))./(repmat(std(data),nr,1));
% compute clusters
clusts=kmeans(data_stnd,5);
% make plots
figure; scatter(data(:,2),data(:,5),50,clusts,'filled');
text(data(:,2),data(:,5),LMEs);
figure; scatter(data(:,1),data(:,6),50,clusts,'filled');
text(data(:,1),data(:,6),LMEs);
figure; scatter(data(:,3),data(:,4),50,clusts,'filled');
text(data(:,3),data(:,4),LMEs);

%% Cluster Analysis of trends
% load and process data
data=readtable('IndsAndStats/TrendsTable-15-Aug-2023.xls');
LMEs = data{:,1};
Vars = data.Properties.VariableNames;
data=table2array(data(:,[2 3 5 6 7 11]));
% Standardize data
nr = size(data,1);
data_stnd = (data - repmat(mean(data),nr,1))./(repmat(std(data),nr,1));
% compute clusters
clusts=kmeans(data_stnd,5);
% make plots
figure; scatter(data(:,2),data(:,5),50,clusts,'filled');
text(data(:,2),data(:,5),LMEs);
figure; scatter(data(:,1),data(:,6),50,clusts,'filled');
text(data(:,1),data(:,6),LMEs);
figure; scatter(data(:,3),data(:,4),50,clusts,'filled');
text(data(:,3),data(:,4),LMEs);

%% Cluster Analysis of amplitudes
% load and process data
data=readtable('IndsAndStats/AmplitudeTable-15-Aug-2023.xls');
LMEs = data{:,1};
Vars = data.Properties.VariableNames;
data=table2array(data(:,[2 3 5 6 7 11]));
% Standardize data
nr = size(data,1);
data_stnd = (data - repmat(mean(data),nr,1))./(repmat(std(data),nr,1));
% compute clusters
clusts=kmeans(data_stnd,5);
% make plots
figure; scatter(data(:,2),data(:,5),50,clusts,'filled');
text(data(:,2),data(:,5),LMEs);
figure; scatter(data(:,1),data(:,6),50,clusts,'filled');
text(data(:,1),data(:,6),LMEs);
figure; scatter(data(:,3),data(:,4),50,clusts,'filled');
text(data(:,3),data(:,4),LMEs);

%% Cluster Analysis of IAVs
% load and process data
data=readtable('IndsAndStats/IAVTable-15-Aug-2023.xls');
LMEs = data{:,1};
Vars = data.Properties.VariableNames;
data=table2array(data(:,[2 3 5 6 7 11]));
% Standardize data
nr = size(data,1);
data_stnd = (data - repmat(mean(data),nr,1))./(repmat(std(data),nr,1));
% compute clusters
clusts=kmeans(data_stnd,5);
% make plots
figure; scatter(data(:,2),data(:,5),50,clusts,'filled');
text(data(:,2),data(:,5),LMEs);
figure; scatter(data(:,1),data(:,6),50,clusts,'filled');
text(data(:,1),data(:,6),LMEs);
figure; scatter(data(:,3),data(:,4),50,clusts,'filled');
text(data(:,3),data(:,4),LMEs);