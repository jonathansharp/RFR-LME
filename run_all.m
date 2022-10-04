% Create surface pCO2 maps in US LMEs using machine learning algorithms
% trained with SOCAT data

% navigate to correct directory
cd /raid/sharp/matlab/US-LMES

% this script grids fCO2 observations from the SOCAT database into grid
% cells of resolution: 0.25 deg lat x 0.25 deg lon x 1 month
grid_vars

% this script loads gridded satellite, model, and reanalysis variables and
% re-grids them to match the size of the fCO2 grids
load_vars

% this script defines predictors variables for algorithm training as X and
% the target variable for algorithm training (i.e. fCO2) as Y
define_x_y

% this script trains self-organizing-maps using the defined predictors to
% cluster data spatiotemporally for algorithm training
cluster_on_grid

% this script trains machine learning algorithms for fCO2 prediction in
% each cluster
fit_algs

% this script applies the machine learning algorithms within each cluster
% to produce fCO2 estimates on the original grids
predict_pco2

% close matlab
exit
