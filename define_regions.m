% Define regions over which to fit fCO2 observations
% 
% Written by J.D. Sharp: 10/21/22
% Last updated by J.D. Sharp: 11/15/22
% 

%% define regions of interest
region = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'HI' ...
          'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};

%% core LMEs

% import LME bounds
lme_shape = shaperead('Data/global/LME66/LMEs66.shp');

% define regions according to numbers
lme_idx.CCS = 10;
lme_idx.GA = 4;
lme_idx.AI = 59;
lme_idx.EBS = 60;
lme_idx.BS = 63;
lme_idx.NBCS = 66;
lme_idx.NE = 12;
lme_idx.SE = 19;
lme_idx.GM = 21;
lme_idx.CS = 25;
lme_idx.HI = 24;

%% Pacific Island LMEs

% import LME grid for US Pacific Islands
lme_grid_is = ncread('Data/global/LME_Scott/mhw_lme_v4.nc','layer');
lme_lon_is = ncread('Data/global/LME_Scott/mhw_lme_v4.nc','longitude');
lme_lat_is = ncread('Data/global/LME_Scott/mhw_lme_v4.nc','latitude');

% define US Pacific Island regions according to numbers
lme_idx_is.AS = 5;
lme_idx_is.JI = 6;
lme_idx_is.PK = 7;
lme_idx_is.HB = 8;
lme_idx_is.JA = 9;
lme_idx_is.WI = 10;
lme_idx_is.GC = 11;

% define LME bounds from grid
for r = 12:18
    idx = bwboundaries(lme_grid_is==lme_idx_is.(region{r}));
    idx = idx{1,1};
    idx_x = idx(:,1);
    idx_y = idx(:,2);
    lme_shape(r+55).X = lme_lon_is(idx_x)';
    lme_shape(r+55).Y = lme_lat_is(idx_y)';
    lme_shape(r+55).LME_NAME = region{r};
    lme_idx.(region{r}) = r+55;
end

clear lme_grid_is lme_lon_is lme_lat_is lme_idx_is r idx_x idx_y idx