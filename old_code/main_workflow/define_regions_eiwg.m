% Define regions over which to fit fCO2 observations
% 
% Written by J.D. Sharp: 5/15/23
% Last updated by J.D. Sharp: 5/15/23
% 

%% define regions of interest
region = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'PI'};

%% Pacific Island LMEs

% import LME bounds
% lme_shape_orig = shaperead('eiwg_boundary_2023/eiwg_boundaries_20230512.shp');

% import LME bounds
lme_shape_tmp = m_shaperead('eiwg_boundary_2023/eiwg_boundaries_20230512');
X = cell(16,1);
Y = cell(16,1);
BoundingBox = cell(16,1);
for r = 1:length(lme_shape_tmp.ncst)
    tmp = lme_shape_tmp.ncst(r);
    tmp = tmp{1};
    X{r} = tmp(:,1)';
    Y{r} = tmp(:,2)';
    BoundingBox{r} = lme_shape_tmp.mbr(r,:);
end
area_km2 = lme_shape_tmp.dbfdata(:,1);
RegionName = lme_shape_tmp.dbfdata(:,2);
lme_shape = struct('BoundingBox',BoundingBox,'X',X,'Y',Y,...
    'area_km2',area_km2,'RegionName',RegionName);
clear lme_shape_tmp tmp X Y BoundingBox area_km2 RegionName

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

% remove weird southern points for NE region
idx = lme_shape(lme_idx.NE).Y < 32;
lme_shape(lme_idx.NE).Y(idx) = [];
lme_shape(lme_idx.NE).X(idx) = [];

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
