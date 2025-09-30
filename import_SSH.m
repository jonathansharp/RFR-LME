% import SSH
function data = import_SSH(dpath,vrs,type,lat,lon,time,yr_end,cmems,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% check for existence of file
if ~isfile(['Data/SSH_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'CMEMS')
    data = import_SSH_CMEMS(dpath,lat,lon,time,yr_end,cmems);
elseif strcmp(type,'ECCO')
    data = import_SSH_ECCO(dpath,lat,lon,time,yr_end);
else
    error('Input variable "type" must be "CMEMS" or "ECCO"');
end

% save data file
ncsave_3d(['Data/SSH_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'SSH' data 'sea surface height anomaly' 'centimeters above something...'});

else

data = ncread(['Data/SSH_' type '_' vrs '.nc'],'SSH');

end

% create ssh animation
if plot_option == 1
    create_animation('SSH',type,time,lat,lon,data,cmocean('haline'),[-1 1],'Sea Surface Height Anomaly','cm');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS SSH
function data = import_SSH_CMEMS(dpath,lat,lon,time,yr_end,cmems)

    % download data if necessary
    rfr_path = pwd; cd(dpath);
    data_id = 'c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1M-m';
    system([cmems.path 'copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id ...
        ' --start-datetime 1998-01-01' ...
        ' --minimum-latitude ' num2str(floor(min(lat))) ...
        ' --maximum-latitude '  num2str(ceil(max(lat))) ...
        ' --minimum-longitude '  num2str(floor(min(lon))) ...
        ' --maximum-longitude '  num2str(ceil(max(lon))) ...
        ' --username ' cmems.usr ' --password ' cmems.pwd]);
    data_id_int = 'c3s_obs-sl_glo_phy-ssh_myint_twosat-l4-duacs-0.25deg_P1M-m';
    system([cmems.path 'copernicusmarine ' ...
        'subset --skip-existing --dataset-id ' data_id_int ...
        ' --end-datetime 2024-12-31 ' ...
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
    data_time = ncread([dpath fpath],'time'); % hours since 1950-01-01
    data_time = datenum(1950,1,1,double(data_time),0,0) + 14; % add 14 days for mid-month

    % read in data
    data = nan(length(data_lon),length(data_lat),(yr_end-1997)*12);
    data_tmp = ncread([dpath fpath],'sla');
    data(:,:,1:length(data_time)) = data_tmp;

    % Interpolate over some gaps in SSH dataset (linear, 1-D, time), then
    % remaining gaps at either end (nearest, 1-D, time)
    for g = 1:length(data_lon)
        for h = 1:length(data_lat)
            if sum(~isnan(squeeze(data(g,h,:)))) >= 200 % check for "too many" NaNs
                % linear interpolation
                idx = ~isnan(squeeze(data(g,h,:)));
                ssh_tmp = interp1(time(idx),squeeze(data(g,h,idx)),time,'linear');
                % then, nearest neighbor interpolation
                idx = ~isnan(ssh_tmp);
                ssh_fit = interp1(time(idx),ssh_tmp(idx),time,'nearest','extrap');
                data(g,h,:) = ssh_fit;
            else
                data(g,h,:) = NaN;
            end
        end
    end

end

% embedded function to import ECCO SSH
function data = import_SSH_ECCO(dpath,lat,lon,time,yr_end)


end

end
