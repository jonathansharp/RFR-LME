% Load North American SOCAT data
% 
% This script loads observations of fCO2 and ancillary variables around
% North America extracted from SOCATv2022 defined by latitude and longitude
% bounds.
% 
% Written by J.D. Sharp: 7/26/22
% Last updated by J.D. Sharp: 1/9/23

%% display status
disp('Downloading SOCAT Data');

%% load regional SOCAT file
% this is obtained by running 'Read_SOCATv3_v2022.m' in the SOCATv2022
% subfolder with the latitude limits of -18 to 82 degrees north and the
% longitude limits of 140 to 302 degrees east
load('data_to_use/SOCATv2022_USA.mat');

%% assemble data into a structure
SOCAT.expocode = Expocode;
SOCAT.latitude = latitude;
SOCAT.longitude = longitude;
SOCAT.temperature = SST;
SOCAT.salinity = sal;
SOCAT.day = day;
SOCAT.month = mon;
SOCAT.year = yr;
% SOCAT.hour = hh;
% SOCAT.minute = mm;
% SOCAT.second = ss;
% SOCAT.pressure = NCEP_SLP;
SOCAT.dist_to_land = dist_to_land;
SOCAT.fCO2 = fCO2rec;
SOCAT.fCO2_flag = fCO2rec_flag;
SOCAT.fCO2_src = fCO2rec_src;
SOCAT.flag = QC_Flag;
% clean up
clearvars -except SOCAT SOCAT_grid region n lme_* 

%% remove observations before 1998
idxyr = SOCAT.year >= 1998;
vars = fieldnames(SOCAT);
for v = 1:numel(vars)
        tempvar = SOCAT.(string(vars(v)));
        SOCAT.(string(vars(v))) = tempvar(idxyr);
end
clear v tempvar idxyr vars

%% determine months since Jan 1 1998
SOCAT.month_since_1998 = ...
    (SOCAT.year-1998).*12 + SOCAT.month;

%% remove observations with flags other than 2 and A/B/C/D
idxflag = SOCAT.fCO2_flag == 2 & ...
    (strcmp(SOCAT.flag,'A') | ...
     strcmp(SOCAT.flag,'B') | ...
     strcmp(SOCAT.flag,'C') | ...
     strcmp(SOCAT.flag,'D'));
vars = fieldnames(SOCAT);
for v = 1:numel(vars)
        tempvar = SOCAT.(char(vars(v)));
        SOCAT.(char(vars(v))) = tempvar(idxflag);
end
clear v idxflag vars tempvar

%% obtain unique integers for each expocode
SOCAT.cruise = ...
    nan(size(SOCAT.expocode));
cruiselist = unique(SOCAT.expocode);
for c = 1:numel(cruiselist)
    idx = strcmp(cruiselist(c),SOCAT.expocode);
    SOCAT.cruise(idx) = c;
end
clear c cruiselist idx

%% visualize number of observations
time = datetime(SOCAT.year,...
    SOCAT.month,SOCAT.day);
figure('visible','on');
histogram(time);
ylabel('Number of observations');
xlabel('Year');
clear time

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/hist_NA.png');
close

%% Save SOCAT structure
if ~isfolder('Data'); mkdir('Data'); end
save('Data/socat_structure','SOCAT','-v7.3');
