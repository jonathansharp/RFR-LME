% import SSH
function data = import_SSH(dpath,vrs,type,lat,lon,time,yr_end,varargin)

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
    data = import_SSH_CMEMS(dpath,lat,lon,time,yr_end);
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
function data = import_SSH_CMEMS(dpath,lat,lon,time,yr_end)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m
    %     > copernicusmarine get --dataset-id cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m
    %     > conda deactivate

    % file path
    fpath = 'CMEMS/SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057/c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1M-m_202411/';
    % fpath = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1M-m_202411/';
    % fpath = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_202112/'; % old
    fpath_int = 'CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.25deg_P1M-m_202311/'; % old

    % load dimensions
    temp_file = '1993/dt_global_twosat_phy_l4_199301_vDT2024-M01.nc';
    % temp_file = '1993/dt_global_allsat_msla_h_y1993_m01.nc'; % old
    inf = ncinfo([dpath fpath temp_file]);
    data_lat = ncread([dpath fpath temp_file],'latitude'); % degrees north
    data_lon = ncread([dpath fpath temp_file],'longitude'); % degrees east
    data_lon_neg = data_lon; data_lon_pos = data_lon;
    data_lon_neg(data_lon_neg >= 0) = NaN;
    data_lon_pos(data_lon_pos < 0) = NaN;
    data_lon_neg = convert_lon(data_lon_neg);
    data_time = ncread([dpath fpath temp_file],'time'); % hours since 1950-01-01
    data_time = datenum(1950,1,1,double(data_time),0,0) + 14; % add 14 days for mid-month

    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-min(lat)));
    [~,idx_maxlat] = min(abs(data_lat-max(lat)));
    [~,idx_minlon1] = min(abs(data_lon_neg-min(lon)));
    [~,idx_maxlon1] = min(abs(data_lon_neg-max(lon)));
    [~,idx_minlon2] = min(abs(data_lon_pos-min(lon)));
    [~,idx_maxlon2] = min(abs(data_lon_pos-max(lon)));
    [~,idx_mintime] = min(abs(data_time-floor(min(time))));
    [~,idx_maxtime] = min(abs(data_time-datenum(yr_end,12,15))); % last month of last year

    % cut down dimensions
    data_lat = data_lat(idx_minlat:idx_maxlat);
    data_lon = [data_lon_pos(idx_minlon2:idx_maxlon2);data_lon_neg(idx_minlon1:idx_maxlon1)];
    data_time = data_time(idx_mintime:idx_maxtime);
   
    % read in data
    date = datevec(time);
    year = date(:,1);
    month = date(:,2);
    for t = 1:(yr_end-1997)*12
        try
            % science-quality data
            data_neglon = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_twosat_phy_l4_' num2str(year(t)) ...
                sprintf('%02d',month(t)) '_vDT2024-M01.nc'],'sla',...
                [idx_minlon1 idx_minlat 1],...
                [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1]);
            data_poslon = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_twosat_phy_l4_' num2str(year(t)) ...
                sprintf('%02d',month(t)) '_vDT2024-M01.nc'],'sla',...
                [idx_minlon2 idx_minlat 1],...
                [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1]);
            data(:,:,t) = cat(1,data_poslon,data_neglon);
            data_time(t) = ncread([dpath fpath num2str(year(t)) ...
                '/dt_global_twosat_phy_l4_' num2str(year(t)) ...
                sprintf('%02d',month(t)) '_vDT2024-M01.nc'],'time');
%             % science-quality data
%             data(:,:,t) = ncread([dpath fpath num2str(year(t)) ...
%                 '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
%                 '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
%                 [idx_minlon idx_minlat 1],...
%                 [1+idx_maxlon-idx_minlon 1+idx_maxlat-idx_minlat 1]);
%             data_time(t) = ncread([dpath fpath num2str(year(t)) ...
%                 '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
%                 '_m' sprintf('%02d',month(t)) '.nc'],'time');
        catch
            try
                % near-real-time data
                data_neg_lon = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                    [idx_minlon1 idx_minlat 1],...
                    [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1]);
                data_pos_lon = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'sla',...
                    [idx_minlon2 idx_minlat 1],...
                    [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1]);
                data(:,:,t) = cat(1,data_poslon,data_neglon);
                data_time(t) = ncread([dpath fpath_int num2str(year(t)) ...
                    '/dt_global_allsat_msla_h_y' num2str(year(t)) ...
                    '_m' sprintf('%02d',month(t)) '.nc'],'time');
            catch
                data(:,:,t) = nan(size(data(:,:,t-1)));
                data_time(t) = data_time(t-1) + 30;
            end
        end
    end

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
