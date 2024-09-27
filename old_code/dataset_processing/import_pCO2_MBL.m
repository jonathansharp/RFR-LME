% Obtain atmospheric pCO2 from NOAA MBL product
function pCO2_atm = import_pCO2_MBL(lat,lon,month,SSS,SST,mslp,ocean_mask,path)

% Open and scan file
file = fopen('data_to_use/MBL_1998_2022.txt');
NOAA_MBL = textscan(file,'%f','Delimiter',',','CommentStyle','#');
fclose(file);

% Reshape according to length of each column
NOAA_MBL = reshape(cell2mat(NOAA_MBL),83,[]);

% Add data to structure
MBL.year = repmat(NOAA_MBL(1,:),(size(NOAA_MBL,1)-1)./2,1);
MBL.CO2 = NOAA_MBL(2:2:end,:);
MBL.err = NOAA_MBL(3:2:end,:);

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
xCO2_atm = nan(length(lat),length(month));

% Interpolate values to latitudes of LME
for m = 1:size(MBL.year_mon,2)
    xCO2_atm(:,m) = interp1(MBL.lat,MBL.CO2_mon(:,m),lat);
end

% % Extend to 2021 using linear interpolation for each month
% for m = 1:12 % for each month
%     for l = 1:size(xCO2_atm,1) % for each latitude
%         xCO2_atm(l,m+size(MBL.year_mon,2)) = ...
%             interp1(MBL.year_mon(1,m:12:size(MBL.year_mon,2)),...
%             xCO2_atm(l,m:12:size(MBL.year_mon,2)),(m-1)/12+2021,'linear','extrap');
%     end
% end

% Replicate across longitudes
xCO2_atm = repmat(permute(xCO2_atm,[3 1 2]),length(lon),1,1);

% Calculate pCO2 from xCO2 (with vapor pressure correction)
% Also inherently eliminates land values
vapor_pressure = vpress(SSS,SST);
pCO2_atm = xCO2_atm.*(mslp-vapor_pressure);

% Remove values outside of ocean mask
for t = 1:length(month)
    pCO2_atm_tmp = pCO2_atm(:,:,t);
    pCO2_atm_tmp(~ocean_mask) = NaN;
    pCO2_atm(:,:,t) = pCO2_atm_tmp;
end

end
