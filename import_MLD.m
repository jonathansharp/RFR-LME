% import MLD
function data_interp = import_MLD(dpath,vrs,type,lat,lon,time,yr_end,cmems,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/MLD_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'CMEMS')
    data_interp = import_MLD_CMEMS(dpath,lat,lon,time,yr_end,cmems);
else
    error('Input variable "type" must be "CMEMS"');
end

% save data file
ncsave_3d(['Data/MLD_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'MLD' data_interp 'mixed layer thickness' 'meters'});

else

data_interp = ncread(['Data/MLD_' type '_' vrs '.nc'],'MLD');

end

% create sst animation
if plot_option == 1
    create_animation('MLD',type,time,lat,lon,data_interp,cmocean('tempo'),[0 200],'Mixed Layer Depth','');
    % create_animation('MLD_anom',type,time,lat,lon,data_interp-mean(data_interp,3,'omitnan'),cmocean('balance'),[-2 2],'Mixed Layer Depth Anomaly','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS MLD
function data_interp = import_MLD_CMEMS(dpath,lat,lon,time,yr_end,cmems)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_myint_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > conda deactivate

    % download data if necessary
    rfr_path = pwd; cd(dpath);
    data_id = 'cmems_mod_glo_phy_my_0.083deg_P1M-m';
    system(['copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id ...
        ' --maximum-depth 1 --start-datetime 1998-01-01' ...
        ' --minimum-latitude ' num2str(floor(min(lat))) ...
        ' --maximum-latitude '  num2str(ceil(max(lat))) ...
        ' --minimum-longitude '  num2str(floor(min(lon))) ...
        ' --maximum-longitude '  num2str(ceil(max(lon))) ...
        ' --username ' cmems.usr ' --password ' cmems.pwd]);
    data_id_int = 'cmems_mod_glo_phy_myint_0.083deg_P1M-m';
    system(['copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id_int ...
        ' --maximum-depth 1 --end-datetime 2024-12-31 ' ...
        ' --minimum-latitude ' num2str(floor(min(lat))) ...
        ' --maximum-latitude '  num2str(ceil(max(lat))) ...
        ' --minimum-longitude '  num2str(floor(min(lon))) ...
        ' --maximum-longitude '  num2str(ceil(max(lon))) ...
        ' --username ' cmems.usr ' --password ' cmems.pwd]);
    cd(rfr_path);

    % file paths
    fpaths = dir(dpath);
    for f = 1:length(fpaths)
        if contains(fpaths(f).name,data_id)
            fpath = fpaths(f).name;
        end
        if contains(fpaths(f).name,data_id_int)
            fpath_int = fpaths(f).name;
        end
    end

    % load dimensions
    % inf = ncinfo([dpath fpath]);
    data_lat = ncread([dpath fpath],'latitude'); % degrees north
    data_lon = ncread([dpath fpath],'longitude'); % degrees east
    data_time_tmp = ncread([dpath fpath],'time'); % hours since 1950-01-01
    data_time_tmp_int = ncread([dpath fpath_int],'time'); % hours since 1950-01-01
    data_time = [data_time_tmp;data_time_tmp_int];
    data_time = datenum(1950,1,1,double(data_time),0,0) + 14; % add 14 days for mid-month

    % read in data
    data = nan(length(data_lon),length(data_lat),(yr_end-1997)*12);
    data_tmp = ncread([dpath fpath],'mlotst');
    data(:,:,1:length(data_time_tmp)) = data_tmp;
    data_tmp_int = ncread([dpath fpath_int],'mlotst');
    data(:,:,length(data_time_tmp)+1:length(data_time_tmp)+...
        length(data_time_tmp_int)) = data_tmp_int;

    % interpolate onto quarter degree grid
    data_interp = nan(length(lon),length(lat),(yr_end-1997)*12);
    [data_lon_grid,data_lat_grid] = ndgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = ndgrid(lon,lat);
    for t = 1:(yr_end-1997)*12
        data_interp(:,:,t) = griddata(double(data_lon_grid),...
            double(data_lat_grid),double(data(:,:,t)),lon_grid,lat_grid);
    end

end

end
