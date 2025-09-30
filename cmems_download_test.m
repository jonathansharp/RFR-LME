% copernicusmarine toolbox download instructions
% https://help.marine.copernicus.eu/en/articles/7970514-copernicus-marine-toolbox-installation
% 
% Do the following to download copernicusmarine and find path:
% 
% conda create --name copernicusmarine conda-forge::copernicusmarine --yes
% conda activate copernicusmarine
% which copernicus marine
% >> ~/.conda/envs/copernicusmarine/bin/copernicusmarine
% conda deactivate

% cmems properties
cmems.path = '/home/sharp/.conda/envs/copernicusmarine/bin/';
cmems.usr = 'jsharp'; cmems.pwd = 'jvqsEZL9'; % fill username and password

% download data if necessary
data_id = 'cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D';
system([cmems.path 'copernicusmarine ' ...
    'subset --dataset-id ' data_id ...
    ' --start-datetime 2025-08-23' ...
    ' --end-datetime 2025-08-23' ...
    ' --minimum-latitude -90' ...
    ' --maximum-latitude 90' ...
    ' --minimum-longitude 0' ...
    ' --maximum-longitude 360' ...
    ' --username ' cmems.usr ' --password ' cmems.pwd]);

% read data
data_file = [data_id '_multi-vars_179.94W-179.94E_89.94S-89.94N_2025-08-23.nc'];
data_lat = ncread(data_file,'latitude'); % degrees north
data_lon = ncread(data_file,'longitude'); % degrees east
data_time = ncread(data_file,'time'); % hours since 1950-01-01
sla = ncread(data_file,'sla'); % sea level altimetry

% plot data
figure;
pcolor(data_lon,data_lat,sla');
title('Sea Level Anomaly');
shading flat; colorbar;
