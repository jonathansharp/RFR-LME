% Define regions over which to fit fCO2 observations
% 
% Written by J.D. Sharp: 5/15/23
% Last updated by J.D. Sharp: 5/15/23
% 

%% define regions of interest
region = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'PI'};
% region = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'HI' ...
%           'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};

%% Pacific Island LMEs

% import LME bounds
lme_shape = shaperead('eiwg_boundary_2023/eiwg_boundaries_20230512.shp');

% define regions according to numbers
lme_idx.CCS = 4;
lme_idx.GA = 7;
lme_idx.AI = 8;
lme_idx.EBS = 10;
lme_idx.BS = 9;
lme_idx.NBCS = 5;
lme_idx.PI = 2;
lme_idx.GM = 3;
lme_idx.CS = 1;
lme_idx.SE = 11;
lme_idx.NE = 6;

% lme_idx.AS = 5;
% lme_idx.JI = 6;
% lme_idx.PK = 7;
% lme_idx.HB = 8;
% lme_idx.JA = 9;
% lme_idx.WI = 10;
% lme_idx.GC = 11;

% %% plot regions for reference
% figure('visible','on');
% worldmap([-18 82],[140 302]);
% box on; hold on;
% for n = 1:length(region)
%     tmp_lon = lme_shape(lme_idx.(region{n})).X';
%     tmp_lat = lme_shape(lme_idx.(region{n})).Y';
%     plotm(tmp_lat,tmp_lon,'k','linewidth',1);
% end
% land = shaperead('usastatehi.shp', 'UseGeoCoords',true);
% geoshow(land,'FaceColor',rgb('grey'));
