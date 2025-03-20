try

% Script to create LME-RFR indicators
vrs = 'SOCATv2024';
dpath = '/raid/Data/';
pred_dims = {'lon' 'lat' 'sin_month_of_year' 'cos_month_of_year' ...
    'year' 'dist' 'bathy'};
pred_vars = {'SSS' 'SSH' 'SST' 'IceC' 'CHL'};
source = {'GLORYS' 'CMEMS' 'OISST' 'OISST' 'CMEMS'};
clust_dims = {};
clust_vars = {'SST' 'SSH'};
num_groups = [3;5;4;4;5;6;4;3;5;3;4];

% this script defines the bounds of the eleven LMEs
[lme_shape,lme_idx,region] = define_lme();
% load_socat(vrs);
% grid_socat(vrs);
import_vars(vrs,dpath,pred_vars,source);
extract_lme(vrs,pred_vars,source,lme_shape,lme_idx,region);
define_x_y(vrs,clust_vars,pred_vars,clust_dims,pred_dims,region);
cluster_lme(vrs,num_groups,region);


catch
    exit
end