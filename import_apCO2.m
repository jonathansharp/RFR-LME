% import atmospheric pCO2
function data_interp = import_apCO2(dpath,vrs,type,lat,lon,time,...
    yr_end,mslp_type,sss_type,sst_type,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

if ~isfile(['Data/apCO2_' type '_' vrs '.nc'])

% Import based on "type"
if strcmp(type,'MBL')
    % import xCO2
    data_interp_temp = import_axCO2_MBL(dpath,lat,lon,time,yr_end);
    % calculate pCO2 from xCO2 (with vapor pressure correction)
    MSLP = ncread(['Data/MSLP_' mslp_type '_' vrs '.nc'],'MSLP');
    SSS = ncread(['Data/SSS_' sss_type '_' vrs '.nc'],'SSS');
    SST = ncread(['Data/SST_' sst_type '_' vrs '.nc'],'SST');
    vapor_pressure = vpress(SSS,SST);
    data_interp = data_interp_temp.*(MSLP-vapor_pressure);
else
    error('Input variable "type" must be "MBL"');
end

% save data file
ncsave_3d(['Data/apCO2_' type '_' vrs '.nc'],{'lon' lon 'longitude' 'degrees east'},...
    {'lat' lat 'latitude' 'degrees north'},...
    {'time' time(1:(yr_end-1997)*12)-datenum(1950,1,1) 'time' 'days since 1950-1-1'},...
    {'apCO2' data_interp 'atmospheric pCO2' 'uatm'});

else

data_interp = ncread(['Data/apCO2_' type '_' vrs '.nc'],'apCO2');

end

% create sst animation
if plot_option == 1
    create_animation('apCO2',type,time,lat,lon,data_interp,cmocean('matter'),[340 420],'Atmospheric pCO2','');
    % create_animation('MSLP_anom',type,time,lat,lon,data_interp-mean(data_interp,3,'omitnan'),cmocean('balance'),[-2 2],'Atmospheric pCO2 Anomaly','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% embedded function to import NOAA MBL xCO2
function data_interp_temp = import_axCO2_MBL(dpath,lat,lon,time,yr_end)
    % file obtained from:
    % 

    % file name
    fnames = dir([dpath 'NOAA-MBL']);
    for n = 1:length(fnames)
        if contains(fnames(n).name,'.txt'); idx = n; end
    end
    fname = [dpath 'NOAA-MBL/' fnames(idx).name];
    
    % Open and scan file
    file = fopen(fname);
    NOAA_MBL = textscan(file,'%f','Delimiter',',','CommentStyle','#');
    fclose(file);
    
    % Reshape according to length of each column
    NOAA_MBL = reshape(cell2mat(NOAA_MBL),83,[]);

    % index based on time
    date = datevec(time);
    [~,idx_mintime] = min(abs(NOAA_MBL(1,:)'-floor(min(date(:,1)))));

    % Add data to structure
    MBL.year = repmat(NOAA_MBL(1,idx_mintime:end),(size(NOAA_MBL,1)-1)./2,1);
    MBL.CO2 = NOAA_MBL(2:2:end,idx_mintime:end);
    MBL.err = NOAA_MBL(3:2:end,idx_mintime:end);
    
    % Define latitudes
    latsin = [-1.00  -0.95  -0.90  -0.85  -0.80  -0.75  -0.70  -0.65  -0.60 ...
              -0.55  -0.50  -0.45  -0.40  -0.35  -0.30  -0.25  -0.20  -0.15 ...
              -0.10  -0.05   0.00   0.05   0.10   0.15   0.20   0.25   0.30 ...
               0.35   0.40   0.45   0.50   0.55   0.60   0.65   0.70   0.75 ...
               0.80   0.85   0.90   0.95   1.00];
    MBL.lat = asind(latsin)';
    
    % Define year fraction interval based on monthly slices
    MBL.year_mon = repmat(min(min(MBL.year)):1/12:max(max(MBL.year)),size(MBL.CO2,1),1);
    MBL.year_mon = MBL.year_mon(:,1:end-1);
    MBL.year_mon = MBL.year_mon + 1/25;
    
    % Interpolate to monthly time slices 
    MBL.CO2_mon = nan(size(MBL.year_mon));
    MBL.err_mon = nan(size(MBL.year_mon));
    for l = 1:size(MBL.CO2,1)
        MBL.CO2_mon(l,:) = interp1(MBL.year(l,:),MBL.CO2(l,:),MBL.year_mon(l,:));
    end
    
    % Pre-allocate
    data_interp_temp = nan(length(lat),(yr_end-1997)*12);
    
    % Interpolate values to latitudes of LME
    for m = 1:size(MBL.year_mon,2)
        data_interp_temp(:,m) = interp1(MBL.lat,MBL.CO2_mon(:,m),lat);
    end
    
    % Extend to last year using linear interpolation for each month if necessary
    if size(MBL.year_mon,2) < (yr_end-1997)*12
        for m = 1:12 % for each month
            for l = 1:size(data_interp_temp,1) % for each latitude
                data_interp_temp(l,m+size(MBL.year_mon,2)) = ...
                    interp1(MBL.year_mon(1,m:12:size(MBL.year_mon,2)),...
                    data_interp_temp(l,m:12:size(MBL.year_mon,2)),(m-1)/12+2023,'linear','extrap');
            end
        end
    end
    
    % Replicate across longitudes
    data_interp_temp = repmat(permute(data_interp_temp,[3 1 2]),length(lon),1,1);

end

end