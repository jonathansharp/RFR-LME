% Run scripts to create surface fCO2 maps in US LMEs using machine learning
% algorithms trained with SOCAT data, create surface TA maps in the same
% LMEs using ESPER algorithms, then calculate ocean acidification
% indicators for each LME.

% add current folder to Matlab path
addpath(genpath(pwd));

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

% define year of data product
rfr_lme_year = 2024;
include_moorings = 1;
moorings_only = 0;

% this script loads SOCAT fCO2 and ancillary data surrounding North
% America extracted from the SOCATv2024 database
% load_socat

% this script grids fCO2 observations from the SOCAT database into grid
% cells of resolution: 0.25 deg lat x 0.25 deg lon x 1 month
% grid_socat

% this script extracts each of eleven LMEs from the gridded data
% surrounding North America
% extract_lme

% this script loads gridded satellite, model, and reanalysis variables and
% re-grids them to match the size of the fCO2 grids
% load_vars

% this script defines predictors variables for algorithm training as X and
% the target variable for algorithm training (i.e. fCO2) as Y
define_x_y

% Set options (determined via 'GMM_test' in 'run_optimization.m')
set_gmm_options

% this script trains Gaussian mixture models using the defined predictors to
% cluster data spatiotemporally for algorithm training
% cluster_on_grid

% Set options ((determined via 'RFR_test' in 'run_optimization.m'))
set_rfr_options

% this script trains machine learning algorithms for fCO2 prediction in
% each cluster
fit_algs_probs

% this script loads error statistics from k-fold algorithm fits and saves
% them in a table
log_errs

% this script applies the machine learning algorithms within each cluster
% to produce fCO2 estimates on the original grids
predict_fCO2_probs

% this script applies the ESPER algorithm to produce TA estimates within
% each LME, then uses fCO2 and TA to calculate OA indicators within each
% LME
% predict_OA

% time series and climatology
% OA_summary_stats
% OA_time_series
% OA_climatology
% region_wide_stats;

% create plots
% plot_full_predictors;
% plot_full_predictors_gif;
% plot_regional_clusters;
% plot_regional_cluster_probabilities;
% plot_regional_deltas;
% plot_full_indicators;
% plot_full_indicators_gif;

% create netCDF files
% matlab_to_netcdf;

% evaluation
% region_wide_stats;
% gridded_product_eval_2;
compare_moorings;
run_evaluation;

% close matlab
exit
