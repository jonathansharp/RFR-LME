%% set parameters
% SOCAT version and data path
vrs = 'SOCATv2025'; dpath = '/raid/Data/';
yr_end = str2num(extractAfter(vrs,'v')) - 1;
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
% probability threshold for model training
thresh = 0.10;

%% set predictor variable sources
source.SSS = 'CMEMS'; source.SSH = 'CMEMS';
source.SST = 'OISST'; source.IceC = 'OISST';
source.CHL = 'NASA'; source.Wind = 'ERA5';
source.MLD = 'CMEMS'; source.MSLP = 'NCEP';
source.apCO2 = 'MBL'; source.Bathy = 'ETOPO';

%% define the bounds of the eleven LMEs
[lme_shape,lme_idx,region] = define_lme();

%% run scripts to create RFR-LME
load_socat(vrs);
grid_socat(vrs,dpath,yr_end);
download_vars(vrs,dpath,source);
import_vars(vrs,dpath,source,yr_end,pred_vars_arc);
extract_lme(vrs,pred_vars,pred_vars_arc,source,lme_shape,lme_idx,region);
define_x_y(vrs,clust_vars,pred_vars,clust_vars_arc,pred_vars_arc,clust_dims,pred_dims,region);
cluster_lme(vrs,num_groups,region,'plot_option',1);
train_rfr(vrs,num_groups,pred_dims,pred_vars,pred_vars_arc,region,thresh,100,2,ceil((2/3)*length(pred_vars)));
predict_fco2(vrs,num_groups,region,'plot_option',1);
calculate_oa(vrs,region,num_groups,'plot_option',1);
matlab_to_netcdf(vrs,lme_shape,lme_idx,region);
regional_stats(vrs,lme_shape,lme_idx,region,yr_end);
create_figures(vrs,lme_shape,lme_idx,region);
log_errs(vrs,num_groups,region);
region_wide_stats(vrs,date);

%% test numbers of clusters
num_groups = repelem((1:10)',11,1);
for ng = 1:size(num_groups,1)
    train_rfr(vrs,num_groups(:,ng),pred_dims,pred_vars,pred_vars_arc,region,...
        thresh,100,2,ceil((2/3)*length(pred_vars)));
end
