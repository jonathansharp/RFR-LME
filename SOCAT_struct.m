% Assemble into structure
SOCAT.(char(region(n))).expocode = Expocode;
SOCAT.(char(region(n))).latitude = latitude;
SOCAT.(char(region(n))).longitude = longitude;
SOCAT.(char(region(n))).temperature = SST;
SOCAT.(char(region(n))).salinity = sal;
SOCAT.(char(region(n))).day = day;
SOCAT.(char(region(n))).month = mon;
SOCAT.(char(region(n))).year = yr;
SOCAT.(char(region(n))).hour = hh;
SOCAT.(char(region(n))).minute = mm;
SOCAT.(char(region(n))).second = ss;
SOCAT.(char(region(n))).pressure = NCEP_SLP;
SOCAT.(char(region(n))).dist_to_land = dist_to_land;
SOCAT.(char(region(n))).fCO2 = fCO2rec;
SOCAT.(char(region(n))).fCO2_flag = fCO2rec_flag;
SOCAT.(char(region(n))).fCO2_src = fCO2rec_src;
SOCAT.(char(region(n))).flag = QC_Flag;

% Clean up
clearvars -except SOCAT region n