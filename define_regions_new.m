% Define regions over which to fit fCO2 observations
% 
% Written by J.D. Sharp: 2/2/23
% Last updated by J.D. Sharp: 2/2/23
% 

%% define regions of interest
region = {'CCS' 'GA' 'EBS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'HI' ...
          'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};

%% Pacific Island LMEs

% import LME grid for US Pacific Islands
lme_grid_is = ncread('Data/global/LME_Mask.nc','layer');
lme_lon_is = ncread('Data/global/LME_Mask.nc','longitude');
lme_lat_is = ncread('Data/global/LME_Mask.nc','latitude');

% define US Pacific Island regions according to numbers
lme_idx_is.GA = 1;
lme_idx_is.EBS = 2;
lme_idx_is.NBCS = 3;
lme_idx_is.HI = 4;
lme_idx_is.AS = 5;
lme_idx_is.JI = 6;
lme_idx_is.PK = 7;
lme_idx_is.HB = 8;
lme_idx_is.JA = 9;
lme_idx_is.WI = 10;
lme_idx_is.GC = 11;
lme_idx_is.CCS = 12;
lme_idx_is.GM = 13;
lme_idx_is.CS = 14;
lme_idx_is.SE = 15;
lme_idx_is.NE = 16;

% define LME bounds from grid
for r = 1:16
    idx = bwboundaries(lme_grid_is==lme_idx_is.(region{r}));
    idx = idx{1,1};
    idx_x = idx(:,1);
    idx_y = idx(:,2);
    lme_shape(r).X = lme_lon_is(idx_x)';
    lme_shape(r).Y = lme_lat_is(idx_y)';
    lme_shape(r).LME_NAME = region{r};
    lme_idx.(region{r}) = r;
end

clear lme_grid_is lme_lon_is lme_lat_is lme_idx_is r idx_x idx_y idx

% plot regions for reference
figure('visible','on');
worldmap([-18 82],[140 302]);
box on; hold on;
for n = 1:length(region)
    tmp_lon = lme_shape(lme_idx.(region{n})).X';
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
land = shaperead('usastatehi.shp', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
