% Create surface pCO2 maps in US LMEs using machine learning algorithms
% trained with SOCAT data

% navigate to correct directory
cd /raid/sharp/matlab/US-LMES

% this script loads SOCAT fCO2 and ancillary data surrounding North
% American extracted from the SOCATv2022 database
% load_socat % ***DONE FOR NOW***

% this script grids fCO2 observations from the SOCAT database into grid
% cells of resolution: 0.25 deg lat x 0.25 deg lon x 1 month
% grid_socat % ***DONE FOR NOW***

% this script extracts each of eleven LMEs from the gridded data
% surrounding North America
% extract_lme % ***DONE FOR NOW***

% this script loads gridded satellite, model, and reanalysis variables and
% re-grids them to match the size of the fCO2 grids
% load_vars % ***DONE FOR NOW***

% this script defines predictors variables for algorithm training as X and
% the target variable for algorithm training (i.e. fCO2) as Y
define_x_y % ***DONE FOR NOW***

% Procedure to test and optimize GMM cluster parameters
% GMM_test
% Set options (determined via GMM_test)
set_gmm_options

% this script trains self-organizing-maps using the defined predictors to
% cluster data spatiotemporally for algorithm training
cluster_on_grid

% Set options (determined via GMM_test)
set_rfr_options

% this script trains machine learning algorithms for fCO2 prediction in
% each cluster
fit_algs_probs

% this script loads error statistics from k-fold algorithm fits and saves
% them in a table
log_errs

% this script applies the machine learning algorithms within each cluster
% to produce fCO2 estimates on the original grids
predict_fco2_probs

% this script applies the ESPER algorithm to produce TA estimates within
% each LME, then uses fCO2 and TA to calculate OA indicators within each
% LME
predict_OA

% close matlab
exit
