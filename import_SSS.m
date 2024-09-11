% import SSS
function data = import_SSS(dpath,type,lat,lon,time)

% Import based on "type"
if strcmp(type,'GLORYS')
    data = import_SSS_GLORYS(type,dpath,lat,lon,time);
elseif strcmp(type,'BASS')
    data = import_SSS_BASS(type,dpath,lat,lon,time);
else
    error('Input variable "type" must be "GLORYS" or "BASS"');
end


% embedded function to import BASS GLORYS
function data = import_SSS_GLORYS(type,dpath,lat,lon,time)
    % files obtained with the copernicusmarine python toolbox:
    %     > cd dpath/CMEMS
    %     > conda activate copernicusmarine
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > copernicusmarine subset --dataset-id cmems_mod_glo_phy_myint_0.083deg_P1M-m --minimum-depth 0 --maximum-depth 0
    %     > conda deactivate

    % file name
    fname = ['cmems_mod_glo_phy_my_0.083deg_P1M-m_multi-vars_180.00W-' ...
        '179.92E_80.00S-90.00N_0.49m_1993-01-01-2021-06-01.nc'];
    fname_int = ['cmems_mod_glo_phy_myint_0.083deg_P1M-m_multi-vars_180.00W-' ...
        '179.92E_80.00S-90.00N_0.49m_2021-07-01-2024-04-01.nc'];
    
    % load dimensions
    inf = ncinfo([dpath 'CMEMS/' fname]);
    data_lat = ncread([dpath 'CMEMS/' fname],'latitude'); % degrees north
    data_lon = ncread([dpath 'CMEMS/' fname],'longitude'); % degrees east
    data_lon_neg = data_lon; data_lon_pos = data_lon;
    data_lon_neg(data_lon_neg >= 0) = NaN;
    data_lon_pos(data_lon_pos < 0) = NaN;
    data_lon_neg = convert_lon(data_lon_neg);
    data_time = ncread([dpath 'CMEMS/' fname],'time'); % hours since 1950-01-01
    data_time = datenum(1950,1,1,data_time,0,0) + 14; % add 14 days for mid-month
    data_time_int = ncread([dpath 'CMEMS/' fname_int],'time'); % hours since 1950-01-01
    data_time_int = datenum(1950,1,1,data_time_int,0,0) + 14; % add 14 days for mid-month
    
    % index based on dimensions
    [~,idx_minlat] = min(abs(data_lat-floor(min(lat))));
    [~,idx_maxlat] = min(abs(data_lat-ceil(max(lat))));
    [~,idx_minlon1] = min(abs(data_lon_neg-floor(min(lon))));
    [~,idx_maxlon1] = min(abs(data_lon_neg-ceil(max(lon))));
    [~,idx_minlon2] = min(abs(data_lon_pos-floor(min(lon))));
    [~,idx_maxlon2] = min(abs(data_lon_pos-ceil(max(lon))));
    [~,idx_mintime] = min(abs(data_time-floor(min(time))));
    [~,idx_maxtime] = min(abs(data_time-ceil(max(time))));
    [~,idx_mintime_int] = min(abs(data_time_int-floor(min(time))));
    [~,idx_maxtime_int] = min(abs(data_time_int-ceil(max(time))));
    
    % cut down dimensions
    data_lat = data_lat(idx_minlat:idx_maxlat);
    data_lon = [data_lon_pos(idx_minlon2:idx_maxlon2);data_lon_neg(idx_minlon1:idx_maxlon1)];
    data_time = data_time(idx_mintime:idx_maxtime);
    data_time_int = data_time_int(idx_mintime_int:idx_maxtime_int);
    
    % read in science-quality data in two longitude chunks
    data_neglon = ncread([dpath 'CMEMS/' fname],'so',[idx_minlon1 idx_minlat 1 idx_mintime],...
        [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1 1+idx_maxtime-idx_mintime]);
    data_poslon = ncread([dpath 'CMEMS/' fname],'so',[idx_minlon2 idx_minlat 1 idx_mintime],...
        [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1 1+idx_maxtime-idx_mintime]);
    
    % read in interim data in two longitude chunks
    data_neglon_int = ncread([dpath 'CMEMS/' fname_int],'so',[idx_minlon1 idx_minlat 1 idx_mintime_int],...
        [1+idx_maxlon1-idx_minlon1 1+idx_maxlat-idx_minlat 1 1+idx_maxtime_int-idx_mintime_int]);
    data_poslon_int = ncread([dpath 'CMEMS/' fname_int],'so',[idx_minlon2 idx_minlat 1 idx_mintime_int],...
        [1+idx_maxlon2-idx_minlon2 1+idx_maxlat-idx_minlat 1 1+idx_maxtime_int-idx_mintime_int]);
    
    % combine data into final grid
    data = cat(4,cat(1,data_poslon,data_neglon),cat(1,data_poslon_int,data_neglon_int));
    data = squeeze(data);
    clear data_poslon data_neglon data_poslon_int data_neglon_int
    
    % interpolate onto quarter degree grid
    data_interp = nan(length(lat),length(lon),length(time));
    [data_lon_grid,data_lat_grid] = meshgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = meshgrid(lon,lat);
    for t = 1:length(time)
        data_interp(:,:,t) = griddata(double(data_lon_grid),...
            double(data_lat_grid),double(data(:,:,t))',lon_grid,lat_grid);
    end

end


% embedded function to import BASS SSS
function data = import_SSS_BASS(type,dpath,lat,lon,time)

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
    data_interp = nan(length(lat),length(lon),length(time));
    [data_lon_grid,data_lat_grid] = meshgrid(data_lon,data_lat);
    [lon_grid,lat_grid] = meshgrid(lon,lat);
    for t = 1:length(time)
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

% create sst animation
create_animation(time,lat,lon,data_interp,cmocean('haline'),[33 37],'Salinity');
% create_animation(time,lat,lon,data_uncer_interp,cmocean('haline'),[33 37],'Salinity');

% save data file
save(['Data/SSS_' type '_' vrs],'data_interp');

% embedded function to create animations
function create_animation(time,lat,lon,z,cmap,colorlims,colorlabel)

    % establish figure
    f = figure; set(f,'color','w','visible','off');

    % loop through months and plot each one
    for m = 1:length(time)
        m_proj('Robinson','lat',[min(lat) max(lat)],'lon',[min(lon) max(lon)]);
        m_pcolor(lon-0.5,lat-0.5,z(:,:,m))
        title(gca,extractAfter(datestr(time(m)),'-'));
        colormap(cmap);
        m_coast('patch',[0.7 0.7 0.7]);
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
        clim(colorlims);
        c=colorbar;
        c.Limits = colorlims;
        c.Label.String = colorlabel;
        c.TickLength = 0;
        % save frame
        if ~isfolder('Figures/SSS/Monthly'); mkdir('Figures/SSS/Monthly'); end
        exportgraphics(f,['Figures/SSS/Monthly/m' num2str(m) '_' strrep(colorlabel,' ','_') '.png']);
        % capture frame
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if m == 1
            imwrite(imind,cm,['Figures/SSS/' strrep(colorlabel,' ','_') '_Animation.gif'],...
                'gif','Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind,cm,['Figures/SSS/' strrep(colorlabel,' ','_') '_Animation.gif'],...
                'gif','WriteMode','append','DelayTime',0.1);
        end
        clf; % clear frame
    end

end

end
