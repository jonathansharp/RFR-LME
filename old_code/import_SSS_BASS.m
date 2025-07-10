% obtain BASS file and NODC climatology file
url = 'https://ftp.cpc.ncep.noaa.gov/precip/BASS/';
fname = 'BASS_V0.Z_MON_1DEG.lnx.B201001';
fname_clim = 'NODC_SAL_1DEG_MON.lnx.CLIM';
websave(['Data/' fname],[url fname]);
websave(['Data/' fname_clim],[url fname_clim]);

% dimensions
nlat = 180;
nlon = 360;
nn = 2; % BASS fields
nd = 24; % climatology depths
lat = linspace(-89.5, 89.5, nlat);
lon = linspace(0.5, 359.5, nlon);

% the size for one month of BASS data:
% 4 bytes (each number) x 360 (lon) x 180 (lat) x 2 (fields, anomaly and error)
fileInfo = dir(['Data/' fname]);
nmonth = fileInfo.bytes / (4* 360 * 180 * 2);
disp(['BASS months = ', num2str(nmonth)]);

% the size for one month of SST climatology data:
% 4 bytes (each number) x 360 (lon) x 180 (lat) x 24 (depth levels)
fileInfo_clim = dir(['Data/' fname_clim]);
nmonth_clim = fileInfo_clim.bytes / (4 * 360 * 180 * 24);
disp(['Climatology months = ', num2str(nmonth_clim)]);

% load sst anomaly data
fid = fopen(['Data/' fname], 'rb');
data = fread(fid, [nlon*nlat*nn nmonth], 'float32');
data = reshape(data,[nlon nlat nn nmonth]);
fclose(fid);
data(data == -999.0) = NaN;
disp(['BASS dimensions: ' num2str(size(data))]);

% load sst climatology data
fid = fopen(['Data/' fname_clim], 'rb');
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

% create sst animation
create_animation(nmonth,lat,lon,squeeze(data(:,:,1,:)),cmocean('haline'),[33 37],'Salinity');
create_animation(nmonth,lat,lon,squeeze(data(:,:,2,:)),cmocean('amp'),[0 1],'Salinity Uncertainty');

% embedded function to create animations
function create_animation(nmonth,lat,lon,z,cmap,colorlims,colorlabel)

    % establish figure
    f = figure; set(f,'color','w','visible','off');

    % loop through months and plot each one
    for m = 1:nmonth
        m_proj('Robinson','lat',[-90 90],'lon',[0 360]);
        m_pcolor(lon-0.5,lat-0.5,z(:,:,m)')
        title(gca,extractAfter(datestr(datenum(2010,m,1)),'-'));
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