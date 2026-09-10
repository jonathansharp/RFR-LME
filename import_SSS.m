% import SSS
function data_interp = import_SSS(dpath,vrs,type,lat,lon,time,yr_end,cmems,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/SSS_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'CMEMS')
    data_interp = import_SSS_CMEMS(dpath,lat,lon,time,yr_end,cmems);
elseif strcmp(type,'BASS')
    data_interp = import_SSS_BASS(dpath,lat,lon,time,yr_end);
else
    error('Input variable "type" must be "CMEMS" or "BASS"');
end

% save data file
ncsave_3d(['Data/SSS_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'SSS' data_interp 'sea surface salinity' ''});

else

data_interp = ncread(['Data/SSS_' type '_' vrs '.nc'],'SSS');

end

% create sst animation
if plot_option == 1
    create_animation('SSS',type,time,lat,lon,data_interp,cmocean('haline'),[32 37],'Salinity','');
    create_animation('SSS_anom',type,time,lat,lon,data_interp-mean(data_interp,3,'omitnan'),cmocean('balance'),[-2 2],'Salinity Anomaly','');
    if strcmp(type,'BASS')
        % create_animation('uSSS',time,lat,lon,data_uncer_interp,cmocean('haline'),[33 37],'Salinity');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS SSS
function data_interp = import_SSS_CMEMS(dpath,lat,lon,time,yr_end,cmems)

     % download data if necessary
    rfr_path = pwd; cd(dpath);
    data_id = 'cmems_mod_glo_phy_my_0.083deg_P1M-m';
    system(['copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id ...
        ' --minimum-depth 0 --maximum-depth 1' ...
        ' --start-datetime 1998-01-01' ...
        ' --minimum-latitude ' num2str(floor(min(lat))) ...
        ' --maximum-latitude '  num2str(ceil(max(lat))) ...
        ' --minimum-longitude '  num2str(floor(min(lon))) ...
        ' --maximum-longitude '  num2str(ceil(max(lon))) ...
        ' --username ' cmems.usr ' --password ' cmems.pwd]);
    data_id_int = 'cmems_mod_glo_phy_myint_0.083deg_P1M-m';
    system(['copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id_int ...
        ' --minimum-depth 0 --maximum-depth 1' ...
        ' --start-datetime 1998-01-01' ...
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
    data_tmp = ncread([dpath fpath],'so');
    data(:,:,1:length(data_time_tmp)) = data_tmp;
    data_tmp_int = ncread([dpath fpath_int],'so');
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


% embedded function to import BASS SSS
function data_interp = import_SSS_BASS(dpath,lat,lon,time,yr_end)

    % obtain BASS file and NODC climatology file
    url = 'https://ftp.cpc.ncep.noaa.gov/precip/BASS/';
    fname = 'BASS_V0.Z_MON_1DEG.lnx.B201001';
    fname_clim = 'NODC_SAL_1DEG_MON.lnx.CLIM';
    websave([dpath 'BASS/' fname],[url fname]);
    websave([dpath 'BASS/' fname_clim],[url fname_clim]);
    
    % dimensions
    nlat = 180;
    nlon = 360;
    nn = 2; % BASS fields
    nd = 24; % climatology depths
    data_lat = linspace(-89.5, 89.5, nlat);
    data_lon = linspace(0.5, 359.5, nlon);
    
    % the size for one month of BASS data:
    % 4 bytes (each number) x 360 (lon) x 180 (lat) x 2 (fields, anomaly and error)
    fileInfo = dir([dpath 'BASS/' fname]);
    nmonth = fileInfo.bytes / (4* 360 * 180 * 2);
    disp(['BASS months = ', num2str(nmonth)]);
    
    % the size for one month of SST climatology data:
    % 4 bytes (each number) x 360 (lon) x 180 (lat) x 24 (depth levels)
    fileInfo_clim = dir([dpath 'BASS/' fname_clim]);
    nmonth_clim = fileInfo_clim.bytes / (4 * 360 * 180 * 24);
    disp(['Climatology months = ', num2str(nmonth_clim)]);
    
    % load sst anomaly data
    fid = fopen([dpath 'BASS/' fname], 'rb');
    data = fread(fid, [nlon*nlat*nn nmonth], 'float32');
    data = reshape(data,[nlon nlat nn nmonth]);
    fclose(fid);
    data(data == -999.0) = NaN;
    disp(['BASS dimensions: ' num2str(size(data))]);
    
    % load sss climatology data
    fid = fopen([dpath 'BASS/' fname_clim], 'rb');
    data_clim = fread(fid, [nlon*nlat*nd nmonth_clim], 'float32');
    data_clim = reshape(data_clim,[nlon nlat nd nmonth_clim]);
    data_clim = squeeze(data_clim(:,:,1,:)); % just the surface layer
    fclose(fid);
    data_clim(data_clim == -999.0) = NaN;
    disp(['Climatology dimensions: ' num2str(size(data_clim))]);
    
    % add BASS to climatology
    data_clim_rep = repmat(data_clim,1,1,ceil(nmonth/12));
    for m = 1:nmonth
        data(:,:,1,m) = data(:,:,1,m) + data_clim_rep(:,:,m);
    end

    % pad data
    year = datevec(time); year = year(:,1);
    max_time = 12*(1+max(year)-2010);
    data = cat(4,nan(nlon,nlat,2,144),data(:,:,:,1:max_time));

    % separate SSS from uncertainty
    data_uncer = squeeze(data(:,:,2,:));
    data = squeeze(data(:,:,1,:));

    % interpolate data onto quarter degree grid
    data_interp = nan(length(lat),length(lon),(yr_end-1997)*12);
    [data_lon_grid,data_lat_grid] = meshgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = meshgrid(lon,lat);
    for t = 1:(yr_end-1997)*12
        data_interp(:,:,t) = griddata(data_lon_grid,...
            data_lat_grid,data(:,:,t)',lon_grid,lat_grid);
    end

    % interpolate uncertainty onto quarter degree grid
    % data_uncer_interp = nan(length(lat),length(lon),length(time));
    % [data_lon_grid,data_lat_grid] = meshgrid(data_lon,data_lat);
    % [lon_grid,lat_grid] = meshgrid(lon,lat);
    % for t = 1:length(time)
    %     data_uncer_interp(:,:,t) = griddata(data_lon_grid,...
    %         data_lat_grid,data_uncer(:,:,t)',lon_grid,lat_grid);
    % end

end

end
