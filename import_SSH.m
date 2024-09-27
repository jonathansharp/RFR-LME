% import SSH
function data = import_SSH(dpath,vrs,type,lat,lon,time,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

% Import based on "type"
if strcmp(type,'CMEMS')
    data = import_SSH_CMEMS(dpath,lat,lon,time);
elseif strcmp(type,'ECCO')
    data = import_SSH_ECCO(dpath,lat,lon,time);
else
    error('Input variable "type" must be "CMEMS" or "ECCO"');
end

% create ssh animation
if plot_option == 1
    create_animation('SSH',type,time,lat,lon,data,cmocean('haline'),[-1 1],'Sea Surface Height Anomaly','cm');
end

% save data file
ncsave_3d(['Data/SSH_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},{'time' time-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'SSH' data 'sea surface height anomaly' 'centimeters above something...'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import CMEMS SSH
function data = import_SSH_CMEMS(dpath,lat,lon,time)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m
    %     > conda deactivate

    % file path
    fpath = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_202112/';
    fpath_int = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m_202311/';

    % load dimensions
    temp_file = '1993/dt_global_allsat_msla_h_y1993_m01.nc';
    inf = ncinfo([dpath fpath temp_file]);
    data_lat = ncread([dpath fpath temp_file],'latitude'); % degrees north
    data_lon = ncread([dpath fpath temp_file],'longitude'); % degrees east
    
    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon] = min(abs(data_lon-min(lon)));
    [~,idx_maxlon] = min(abs(data_lon-max(lon)));

    % cut down dimensions
    data_lat = data_lat(idx_minlat:idx_maxlat);
    data_lon = data_lon(idx_minlon:idx_maxlon);

    % read in data
    date = datevec(time);
    year = date(:,1);
    month = date(:,2);
    for t = 1:length(time)
        try
            % science-quality data
            data(:,:,t) = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                [idx_minlon idx_minlat 1],...
                [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1]);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                '_m' sprintf('%02d',month(t)) '.nc'],'time');
        catch
            try
                % near-real-time data
                data(:,:,t) = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                    [idx_minlon idx_minlat 1],...
                    [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1]);
                data_time(t) = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'time');
            catch
            end
        end
    end

end

% embedded function to import ECCO SSH
function data_interp = import_SSH_ECCO(dpath,lat,lon,time)


end

end
