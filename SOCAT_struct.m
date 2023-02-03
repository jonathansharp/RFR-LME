% Assemble into structure
SOCAT.expocode = Expocode;
SOCAT.latitude = latitude;
SOCAT.longitude = longitude;
SOCAT.temperature = SST;
SOCAT.salinity = sal;
SOCAT.day = day;
SOCAT.month = mon;
SOCAT.year = yr;
SOCAT.hour = hh;
SOCAT.minute = mm;
SOCAT.second = ss;
SOCAT.pressure = NCEP_SLP;
SOCAT.dist_to_land = dist_to_land;
SOCAT.fCO2 = fCO2rec;
SOCAT.fCO2_flag = fCO2rec_flag;
SOCAT.fCO2_src = fCO2rec_src;
SOCAT.flag = QC_Flag;

% Clean up
clearvars -except SOCAT SOCAT_grid region n lme_shape lme_idx