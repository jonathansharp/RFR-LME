%% set user-defined parameters
% SOCAT version, data paths, and passwords
paths = rfr_lme_paths(); % this helper set paths to data and figures, edit the helper function to change paths
vrs = 'SOCATv2026'; % change this version tag annually
cmems.path = ''; % add path to where copernicusmarine toolbox is located
cmems.usr = ''; % add username for copernicus marine data system
cmems.pwd = ''; % add password for copernicus marine data system
nasa.token = ''; % add personal token for NASA Earthdata
cds.token =''; % add personal token for Climate Data Store

%% YOU SHOULDN'T NEED TO EDIT ANYTHING BELOW HERE TO RUN CODE %%

%% set additional parameters
% Coordinates and variables to be used for models
pred_dims = {'lon' 'lat' 'sin_month_of_year' 'cos_month_of_year' 'year' 'dist'};
pred_vars = {'SSS' 'SSH' 'SST' 'IceC' 'CHL' 'Wind' 'MLD' 'MSLP' 'apCO2' 'Bathy'};
pred_vars_arc = {'SSS' 'SST' 'IceC' 'Wind' 'MLD' 'MSLP' 'apCO2' 'Bathy'};
% Coordinates and variables to be used for clustering
clust_dims = {};
clust_vars = {'SST' 'MSLP' 'CHL'};
clust_vars_arc = {'SST' 'MSLP' 'Wind'};
% number of rfr groups
num_groups = [3;5;4;4;5;6;4;3;5;3;4];
yr_end = str2double(extractAfter(vrs,'v')) - 1;
% probability threshold for model training
thresh = 0.10; % probability threshold used to determine whether a data point belongs in a given cluster
clust_opt = 'var'; % set to 'var' to cluster based on variability over time (static clusters) and 'abs' to cluster based on absolute values (time-varying clusters)

%% set predictor variable sources
source.SSS = 'CMEMS'; source.SSH = 'CMEMS';
source.SST = 'OISST'; source.IceC = 'OISST';
source.CHL = 'NASA'; source.Wind = 'ERA5';
source.MLD = 'CMEMS'; source.MSLP = 'NCEP';
source.apCO2 = 'MBL'; source.Bathy = 'ETOPO';

%% define the bounds of the eleven LMEs
[lme_shape,lme_idx,region] = define_lme(paths);

%% run scripts to create RFR-LME
load_socat(vrs,paths);
grid_socat(vrs,paths.dpath,yr_end,paths);
import_vars(vrs,paths.dpath,source,yr_end,pred_vars_arc,cmems,nasa,cds,paths);
extract_lme(vrs,pred_vars,pred_vars_arc,source,lme_shape,lme_idx,region,paths);
define_x_y(vrs,clust_vars,pred_vars,clust_vars_arc,pred_vars_arc,clust_dims,pred_dims,region,paths);
cluster_lme(vrs,num_groups,region,paths,'plot_option',0,'cluster_option',clust_opt);
train_rfr(vrs,num_groups,pred_dims,pred_vars,pred_vars_arc,region,...
    thresh,100,2,ceil((2/3)*length(pred_vars)),paths);
log_errs(vrs,num_groups,region,paths);
predict_fco2(vrs,num_groups,region,paths,'plot_option',1);
calculate_oa(vrs,region,num_groups,paths,'plot_option',1);
matlab_to_netcdf_new(vrs,region,paths);
regional_stats(vrs,lme_shape,lme_idx,region,yr_end,paths);
create_figures(vrs,lme_shape,lme_idx,region,paths);
region_wide_stats(vrs,date,paths);

%% test numbers of clusters
num_groups = repelem((1:10),11,1);
for ng = 1:size(num_groups,1)
    train_rfr(vrs,num_groups(:,ng),pred_dims,pred_vars,pred_vars_arc,region,...
        thresh,100,2,ceil((2/3)*length(pred_vars)),paths);
end
