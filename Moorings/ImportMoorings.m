%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import mooring data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set paths and pre-allocate variables
% filePath = matlab.desktop.editor.getActiveFilename;
% filePath = erase(filePath,'ImportMoorings.m');
filePath = '/raid/sharp/matlab/US-LMES/Moorings/';
fileinfo = dir(strcat(filePath,'*.csv'));
fnames = {fileinfo.name}';
structnames = extractBefore(fnames,'_');
fnames2 = strrep(fnames,'-','_');
structnames = strrep(structnames,'-','_');
warning('OFF','MATLAB:table:ModifiedAndSavedVarnames');
name = []; lat = []; lon = []; date = [];
time = []; sst = []; sal = [];
fCO2SW = []; fCO2Air = [];

%% CCE1
MOORING.CCE1.name    = name;
MOORING.CCE1.lat     = lat;
MOORING.CCE1.lon     = lon;
MOORING.CCE1.date    = date;
MOORING.CCE1.time    = time;
MOORING.CCE1.sst     = sst;
MOORING.CCE1.sal     = sal;
MOORING.CCE1.fCO2SW  = fCO2SW;
MOORING.CCE1.fCO2Air = fCO2Air;
idx = find(strcmp('CCE1',structnames));
for n = 1:size(idx,1)
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'VariableNamingRule','modify');
    if n == 11
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.CCE1.name    = [MOORING.CCE1.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.CCE1.lat     = [MOORING.CCE1.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.CCE1.lon     = [MOORING.CCE1.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.CCE1.date    = [MOORING.CCE1.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.CCE1.time   = [MOORING.CCE1.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.CCE1.time   = [MOORING.CCE1.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.CCE1.sst     = [MOORING.CCE1.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.CCE1.sal     = [MOORING.CCE1.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    if n == 4 || n == 7 || n == 8
        MOORING.CCE1.fCO2SW  = [MOORING.CCE1.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat__uatm_(idxQC)];
        MOORING.CCE1.fCO2Air = [MOORING.CCE1.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat__uatm_(idxQC)];
    else
        MOORING.CCE1.fCO2SW  = [MOORING.CCE1.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
        MOORING.CCE1.fCO2Air = [MOORING.CCE1.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
    end
end
MOORING.CCE1.date(MOORING.CCE1.date(:,1)<1000)=MOORING.CCE1.date(MOORING.CCE1.date(:,1)<1000)+2000;
MOORING.CCE1.sst(MOORING.CCE1.sst==-999) = NaN;
MOORING.CCE1.sal(MOORING.CCE1.sal==-999) = NaN;
MOORING.CCE1.fCO2SW(MOORING.CCE1.fCO2SW==-999) = NaN;
MOORING.CCE1.fCO2Air(MOORING.CCE1.fCO2Air==-999) = NaN;
MOORING.CCE1.lon(MOORING.CCE1.lon<0) = MOORING.CCE1.lon(MOORING.CCE1.lon<0) + 360;

%% CCE2
MOORING.CCE2.name = name;
MOORING.CCE2.lat  = lat;
MOORING.CCE2.lon  = lon;
MOORING.CCE2.date = date;
MOORING.CCE2.time = time;
MOORING.CCE2.sst  = sst;
MOORING.CCE2.sal  = sal;
MOORING.CCE2.fCO2SW  = fCO2SW;
MOORING.CCE2.fCO2Air = fCO2Air;
idx = find(strcmp('CCE2',structnames));
for n = 1:size(idx,1)
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',4,'VariableNamingRule','modify');
    if n == 1
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QFXCO2SW_wet_ == 2;
    else
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.CCE2.name    = [MOORING.CCE2.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.CCE2.lat     = [MOORING.CCE2.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.CCE2.lon     = [MOORING.CCE2.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.CCE2.date    = [MOORING.CCE2.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.CCE2.time   = [MOORING.CCE2.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.CCE2.time   = [MOORING.CCE2.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.CCE2.sst     = [MOORING.CCE2.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.CCE2.sal     = [MOORING.CCE2.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    if n == 3 || n == 8 || n == 11 || n == 12
        MOORING.CCE2.fCO2SW  = [MOORING.CCE2.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat__uatm_(idxQC)];
        MOORING.CCE2.fCO2Air = [MOORING.CCE2.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat__uatm_(idxQC)];
    else
        MOORING.CCE2.fCO2SW  = [MOORING.CCE2.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
        MOORING.CCE2.fCO2Air = [MOORING.CCE2.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
    end
end
MOORING.CCE2.date(MOORING.CCE2.date(:,1)<1000)=MOORING.CCE2.date(MOORING.CCE2.date(:,1)<1000)+2000;
MOORING.CCE2.sst(MOORING.CCE2.sst==-999) = NaN;
MOORING.CCE2.sal(MOORING.CCE2.sal==-999) = NaN;
MOORING.CCE2.fCO2SW(MOORING.CCE2.fCO2SW==-999) = NaN;
MOORING.CCE2.fCO2Air(MOORING.CCE2.fCO2Air==-999) = NaN;
MOORING.CCE2.lon(MOORING.CCE2.lon<0) = MOORING.CCE2.lon(MOORING.CCE2.lon<0) + 360;

%% Southeast Alaska
MOORING.Southeast.name = name;
MOORING.Southeast.lat  = lat;
MOORING.Southeast.lon  = lon;
MOORING.Southeast.date = date;
MOORING.Southeast.time = time;
MOORING.Southeast.sst  = sst;
MOORING.Southeast.sal  = sal;
MOORING.Southeast.fCO2SW  = fCO2SW;
MOORING.Southeast.fCO2Air = fCO2Air;
idx = find(strcmp('Southeast',structnames));
for n = 1:size(idx,1)
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'VariableNamingRule','modify');
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    MOORING.Southeast.name    = [MOORING.Southeast.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.Southeast.lat     = [MOORING.Southeast.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.Southeast.lon     = [MOORING.Southeast.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.Southeast.date    = [MOORING.Southeast.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.Southeast.time   = [MOORING.Southeast.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.Southeast.time   = [MOORING.Southeast.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.Southeast.sst     = [MOORING.Southeast.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.Southeast.sal     = [MOORING.Southeast.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.Southeast.fCO2SW  = [MOORING.Southeast.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.Southeast.fCO2Air = [MOORING.Southeast.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.Southeast.date(MOORING.Southeast.date(:,1)<1000)=MOORING.Southeast.date(MOORING.Southeast.date(:,1)<1000)+2000;
MOORING.Southeast.sst(MOORING.Southeast.sst==-999) = NaN;
MOORING.Southeast.sal(MOORING.Southeast.sal==-999) = NaN;
MOORING.Southeast.fCO2SW(MOORING.Southeast.fCO2SW==-999) = NaN;
MOORING.Southeast.fCO2Air(MOORING.Southeast.fCO2Air==-999) = NaN;
MOORING.Southeast.lon(MOORING.Southeast.lon<0) = MOORING.Southeast.lon(MOORING.Southeast.lon<0) + 360;

%% NH10
MOORING.NH10.name = name;
MOORING.NH10.lat  = lat;
MOORING.NH10.lon  = lon;
MOORING.NH10.date = date;
MOORING.NH10.time = time;
MOORING.NH10.sst  = sst;
MOORING.NH10.sal  = sal;
MOORING.NH10.fCO2SW  = fCO2SW;
MOORING.NH10.fCO2Air = fCO2Air;
idx = find(strcmp('NH10',structnames));
for n = 1:size(idx,1)
    MOORING.(extractBefore(fnames{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames2{idx(n)}),'VariableNamingRule','modify');
    if n == 1
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.NH10.name    = [MOORING.NH10.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.NH10.lat     = [MOORING.NH10.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.NH10.lon     = [MOORING.NH10.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.NH10.date    = [MOORING.NH10.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.NH10.time   = [MOORING.NH10.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.NH10.time   = [MOORING.NH10.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.NH10.sst     = [MOORING.NH10.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.NH10.sal     = [MOORING.NH10.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.NH10.fCO2SW  = [MOORING.NH10.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.NH10.fCO2Air = [MOORING.NH10.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.NH10.date(MOORING.NH10.date(:,1)<1000)=MOORING.NH10.date(MOORING.NH10.date(:,1)<1000)+2000;
MOORING.NH10.sst(MOORING.NH10.sst==-999) = NaN;
MOORING.NH10.sal(MOORING.NH10.sal==-999) = NaN;
MOORING.NH10.fCO2SW(MOORING.NH10.fCO2SW==-999) = NaN;
MOORING.NH10.fCO2Air(MOORING.NH10.fCO2Air==-999) = NaN;
MOORING.NH10.lon(MOORING.NH10.lon<0) = MOORING.NH10.lon(MOORING.NH10.lon<0) + 360;

%% Cape Arrago
MOORING.CB_06.name = name;
MOORING.CB_06.lat  = lat;
MOORING.CB_06.lon  = lon;
MOORING.CB_06.date = date;
MOORING.CB_06.time = time;
MOORING.CB_06.sst  = sst;
MOORING.CB_06.sal  = sal;
MOORING.CB_06.fCO2SW  = fCO2SW;
MOORING.CB_06.fCO2Air = fCO2Air;
idx = find(strcmp('CB_06',structnames));
for n = 1:size(idx,1)
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'VariableNamingRule','modify');
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    MOORING.CB_06.name    = [MOORING.CB_06.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.CB_06.lat     = [MOORING.CB_06.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.CB_06.lon     = [MOORING.CB_06.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.CB_06.date    = [MOORING.CB_06.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.CB_06.time   = [MOORING.CB_06.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.CB_06.time   = [MOORING.CB_06.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.CB_06.sst     = [MOORING.CB_06.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.CB_06.sal     = [MOORING.CB_06.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    if n == 2 || n == 3 || n == 5 || n == 6
        MOORING.CB_06.fCO2SW  = [MOORING.CB_06.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat__uatm_(idxQC)];
        MOORING.CB_06.fCO2Air = [MOORING.CB_06.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat__uatm_(idxQC)];
    else
        MOORING.CB_06.fCO2SW  = [MOORING.CB_06.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
        MOORING.CB_06.fCO2Air = [MOORING.CB_06.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
    end
end
MOORING.CB_06.date(MOORING.CB_06.date(:,1)<1000)=MOORING.CB_06.date(MOORING.CB_06.date(:,1)<1000)+2000;
MOORING.CB_06.sst(MOORING.CB_06.sst==-999) = NaN;
MOORING.CB_06.sal(MOORING.CB_06.sal==-999) = NaN;
MOORING.CB_06.fCO2SW(MOORING.CB_06.fCO2SW==-999) = NaN;
MOORING.CB_06.fCO2Air(MOORING.CB_06.fCO2Air==-999) = NaN;
MOORING.CB_06.lon(MOORING.CB_06.lon<0) = MOORING.CB_06.lon(MOORING.CB_06.lon<0) + 360;

%% Cape Elizabeth
MOORING.WA.name = name;
MOORING.WA.lat  = lat;
MOORING.WA.lon  = lon;
MOORING.WA.date = date;
MOORING.WA.time = time;
MOORING.WA.sst  = sst;
MOORING.WA.sal  = sal;
MOORING.WA.fCO2SW  = fCO2SW;
MOORING.WA.fCO2Air = fCO2Air;
idx = find(strcmp('WA',structnames));
for n = 1:size(idx,1)
    HL = 4;
    if n == 6 || n == 7 || n == 10 || n == 11 || n == 13 || n == 14
        HL = 0;
    end
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
    if n == 6 || n == 7 || n == 10 || n == 11 || n == 13 || n == 14
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.WA.name    = [MOORING.WA.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.WA.lat     = [MOORING.WA.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.WA.lon     = [MOORING.WA.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.WA.date    = [MOORING.WA.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.WA.time   = [MOORING.WA.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.WA.time   = [MOORING.WA.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.WA.sst     = [MOORING.WA.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.WA.sal     = [MOORING.WA.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.WA.fCO2SW  = [MOORING.WA.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.WA.fCO2Air = [MOORING.WA.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.WA.date(MOORING.WA.date(:,1)<1000)=MOORING.WA.date(MOORING.WA.date(:,1)<1000)+2000;
MOORING.WA.sst(MOORING.WA.sst==-999) = NaN;
MOORING.WA.sal(MOORING.WA.sal==-999) = NaN;
MOORING.WA.fCO2SW(MOORING.WA.fCO2SW==-999) = NaN;
MOORING.WA.fCO2Air(MOORING.WA.fCO2Air==-999) = NaN;
MOORING.WA.lon(MOORING.WA.lon<0) = MOORING.WA.lon(MOORING.WA.lon<0) + 360;

%% Cha ba
MOORING.LaPush.name = name;
MOORING.LaPush.lat  = lat;
MOORING.LaPush.lon  = lon;
MOORING.LaPush.date = date;
MOORING.LaPush.time = time;
MOORING.LaPush.sst  = sst;
MOORING.LaPush.sal  = sal;
MOORING.LaPush.fCO2SW  = fCO2SW;
MOORING.LaPush.fCO2Air = fCO2Air;
idx = find(strcmp('LaPush',structnames));
for n = 1:size(idx,1)
    HL = 4; if n == 2; HL = 0; end
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
    if n == 7
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF_xCO2_SW == 2;
    else
    idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.LaPush.name    = [MOORING.LaPush.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.LaPush.lat     = [MOORING.LaPush.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.LaPush.lon     = [MOORING.LaPush.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.LaPush.date    = [MOORING.LaPush.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.LaPush.time   = [MOORING.LaPush.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.LaPush.time   = [MOORING.LaPush.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.LaPush.sst     = [MOORING.LaPush.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.LaPush.sal     = [MOORING.LaPush.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    if n == 3 || n == 5 || n == 12 || n == 15
        MOORING.LaPush.fCO2SW  = [MOORING.LaPush.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat__uatm_(idxQC)];
        MOORING.LaPush.fCO2Air = [MOORING.LaPush.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat__uatm_(idxQC)];
    else
        MOORING.LaPush.fCO2SW  = [MOORING.LaPush.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
        MOORING.LaPush.fCO2Air = [MOORING.LaPush.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
    end
end
MOORING.LaPush.date(MOORING.LaPush.date(:,1)<1000)=MOORING.LaPush.date(MOORING.LaPush.date(:,1)<1000)+2000;
MOORING.LaPush.sst(MOORING.LaPush.sst==-999) = NaN;
MOORING.LaPush.sal(MOORING.LaPush.sal==-999) = NaN;
MOORING.LaPush.fCO2SW(MOORING.LaPush.fCO2SW==-999) = NaN;
MOORING.LaPush.fCO2Air(MOORING.LaPush.fCO2Air==-999) = NaN;
MOORING.LaPush.lon(MOORING.LaPush.lon<0) = MOORING.LaPush.lon(MOORING.LaPush.lon<0) + 360;

%% Kodiak
MOORING.Kodiak.name = name;
MOORING.Kodiak.lat  = lat;
MOORING.Kodiak.lon  = lon;
MOORING.Kodiak.date = date;
MOORING.Kodiak.time = time;
MOORING.Kodiak.sst  = sst;
MOORING.Kodiak.sal  = sal;
MOORING.Kodiak.fCO2SW  = fCO2SW;
MOORING.Kodiak.fCO2Air = fCO2Air;
idx = find(strcmp('Kodiak',structnames));
for n = 1:size(idx,1)
    HL = 4;
    if n == 2 || n == 3 || n == 4
        HL = 0;
    end
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
    if n == 1 || n == 2 || n == 3 || n == 4
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.Kodiak.name    = [MOORING.Kodiak.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.Kodiak.lat     = [MOORING.Kodiak.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.Kodiak.lon     = [MOORING.Kodiak.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.Kodiak.date    = [MOORING.Kodiak.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.Kodiak.time   = [MOORING.Kodiak.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.Kodiak.time   = [MOORING.Kodiak.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.Kodiak.sst     = [MOORING.Kodiak.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.Kodiak.sal     = [MOORING.Kodiak.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.Kodiak.fCO2SW  = [MOORING.Kodiak.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.Kodiak.fCO2Air = [MOORING.Kodiak.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.Kodiak.date(MOORING.Kodiak.date(:,1)<1000)=MOORING.Kodiak.date(MOORING.Kodiak.date(:,1)<1000)+2000;
MOORING.Kodiak.sst(MOORING.Kodiak.sst==-999) = NaN;
MOORING.Kodiak.sal(MOORING.Kodiak.sal==-999) = NaN;
MOORING.Kodiak.fCO2SW(MOORING.Kodiak.fCO2SW==-999) = NaN;
MOORING.Kodiak.fCO2Air(MOORING.Kodiak.fCO2Air==-999) = NaN;
MOORING.Kodiak.lon(MOORING.Kodiak.lon<0) = MOORING.Kodiak.lon(MOORING.Kodiak.lon<0) + 360;

%% GAKOA
MOORING.GAKOA.name = name;
MOORING.GAKOA.lat  = lat;
MOORING.GAKOA.lon  = lon;
MOORING.GAKOA.date = date;
MOORING.GAKOA.time = time;
MOORING.GAKOA.sst  = sst;
MOORING.GAKOA.sal  = sal;
MOORING.GAKOA.fCO2SW  = fCO2SW;
MOORING.GAKOA.fCO2Air = fCO2Air;
idx = find(strcmp('GAKOA',structnames));
for n = 1:size(idx,1)
    HL = 4;
    if n == 3 || n == 4 || n == 5 || n == 10 || n == 11 || n == 12
        HL = 0;
    end
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
    if n == 3 || n == 4 || n == 5 || n == 10 || n == 11 || n == 12
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.GAKOA.name    = [MOORING.GAKOA.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.GAKOA.lat     = [MOORING.GAKOA.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.GAKOA.lon     = [MOORING.GAKOA.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.GAKOA.date    = [MOORING.GAKOA.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.GAKOA.time   = [MOORING.GAKOA.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.GAKOA.time   = [MOORING.GAKOA.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.GAKOA.sst     = [MOORING.GAKOA.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.GAKOA.sal     = [MOORING.GAKOA.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.GAKOA.fCO2SW  = [MOORING.GAKOA.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.GAKOA.fCO2Air = [MOORING.GAKOA.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.GAKOA.date(MOORING.GAKOA.date(:,1)<1000)=MOORING.GAKOA.date(MOORING.GAKOA.date(:,1)<1000)+2000;
MOORING.GAKOA.sst(MOORING.GAKOA.sst==-999) = NaN;
MOORING.GAKOA.sal(MOORING.GAKOA.sal==-999) = NaN;
MOORING.GAKOA.fCO2SW(MOORING.GAKOA.fCO2SW==-999) = NaN;
MOORING.GAKOA.fCO2Air(MOORING.GAKOA.fCO2Air==-999) = NaN;
MOORING.GAKOA.lon(MOORING.GAKOA.lon<0) = MOORING.GAKOA.lon(MOORING.GAKOA.lon<0) + 360;

%% M2
MOORING.M2.name = name;
MOORING.M2.lat  = lat;
MOORING.M2.lon  = lon;
MOORING.M2.date = date;
MOORING.M2.time = time;
MOORING.M2.sst  = sst;
MOORING.M2.sal  = sal;
MOORING.M2.fCO2SW  = fCO2SW;
MOORING.M2.fCO2Air = fCO2Air;
idx = find(strcmp('M2',structnames));
for n = 1:size(idx,1)
    HL = 4;
    if n == 3 || n == 4
        HL = 0;
    end
    MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
        readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
    if n == 3 || n == 4
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
    else
        idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
    end
    MOORING.M2.name    = [MOORING.M2.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
    MOORING.M2.lat     = [MOORING.M2.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
    MOORING.M2.lon     = [MOORING.M2.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
    MOORING.M2.date    = [MOORING.M2.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
    if ~iscell(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC))
        MOORING.M2.time   = [MOORING.M2.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
    else
        hr = cellfun(@str2num,extractBefore(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        mn = cellfun(@str2num,extractAfter(MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC),':'));
        MOORING.M2.time   = [MOORING.M2.time;duration(hr,mn,zeros(size(hr)))];
    end
    MOORING.M2.sst     = [MOORING.M2.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
    MOORING.M2.sal     = [MOORING.M2.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
    MOORING.M2.fCO2SW  = [MOORING.M2.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
    MOORING.M2.fCO2Air = [MOORING.M2.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
end
MOORING.M2.date(MOORING.M2.date(:,1)<1000)=MOORING.M2.date(MOORING.M2.date(:,1)<1000)+2000;
MOORING.M2.sst(MOORING.M2.sst==-999) = NaN;
MOORING.M2.sal(MOORING.M2.sal==-999) = NaN;
MOORING.M2.fCO2SW(MOORING.M2.fCO2SW==-999) = NaN;
MOORING.M2.fCO2Air(MOORING.M2.fCO2Air==-999) = NaN;
MOORING.M2.lon(MOORING.M2.lon<0) = MOORING.M2.lon(MOORING.M2.lon<0) + 360;

%% Kwakshua Channel 
% MOORING.KwakshuaChannel.name = name;
% MOORING.KwakshuaChannel.lat  = lat;
% MOORING.KwakshuaChannel.lon  = lon;
% MOORING.KwakshuaChannel.date = date;
% MOORING.KwakshuaChannel.time = time;
% MOORING.KwakshuaChannel.sst  = sst;
% MOORING.KwakshuaChannel.sal  = sal;
% MOORING.KwakshuaChannel.fCO2SW  = fCO2SW;
% MOORING.KwakshuaChannel.fCO2Air = fCO2Air;
% idx = find(strcmp('187F20180501',structnames));
% for n = 1:size(idx,1)
%     HL = 6;
%     MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))) = ...
%         readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
%     idxQC = MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).CO2SWQF == 2;
%     MOORING.KwakshuaChannel.name    = [MOORING.KwakshuaChannel.name;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).MooringName(idxQC)];
%     MOORING.KwakshuaChannel.lat     = [MOORING.KwakshuaChannel.lat;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).Latitude(idxQC)];
%     MOORING.KwakshuaChannel.lon     = [MOORING.KwakshuaChannel.lon;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).Longitude(idxQC)];
%     MOORING.KwakshuaChannel.date    = [MOORING.KwakshuaChannel.date;datevec(MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).Date(idxQC))];
%     if ~iscell(MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).TimeUTC(idxQC))
%         MOORING.KwakshuaChannel.time   = [MOORING.KwakshuaChannel.time;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).TimeUTC(idxQC)];
%     else
%         hr = cellfun(@str2num,extractBefore(MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).TimeUTC(idxQC),':'));
%         mn = cellfun(@str2num,extractAfter(MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).TimeUTC(idxQC),':'));
%         MOORING.KwakshuaChannel.time   = [MOORING.KwakshuaChannel.time;duration(hr,mn,zeros(size(hr)))];
%     end
%     MOORING.KwakshuaChannel.sst     = [MOORING.KwakshuaChannel.sst;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).SST_C_(idxQC)];
%     MOORING.KwakshuaChannel.sal     = [MOORING.KwakshuaChannel.sal;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).Salinity(idxQC)];
%     MOORING.KwakshuaChannel.fCO2SW  = [MOORING.KwakshuaChannel.fCO2SW;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).fCO2SW_sat_Uatm(idxQC)];
%     MOORING.KwakshuaChannel.fCO2Air = [MOORING.KwakshuaChannel.fCO2Air;MOORING.(strcat('b',extractBefore(fnames2{idx(n)},'.'))).fCO2Air_sat_Uatm(idxQC)];
% end
% MOORING.KwakshuaChannel.date(MOORING.KwakshuaChannel.date(:,1)<1000)=MOORING.KwakshuaChannel.date(MOORING.KwakshuaChannel.date(:,1)<1000)+2000;
% MOORING.KwakshuaChannel.sst(MOORING.KwakshuaChannel.sst==-999) = NaN;
% MOORING.KwakshuaChannel.sal(MOORING.KwakshuaChannel.sal==-999) = NaN;
% MOORING.KwakshuaChannel.fCO2SW(MOORING.KwakshuaChannel.fCO2SW==-999) = NaN;
% MOORING.KwakshuaChannel.fCO2Air(MOORING.KwakshuaChannel.fCO2Air==-999) = NaN;
% MOORING.KwakshuaChannel.lon(MOORING.KwakshuaChannel.lon<0) = MOORING.KwakshuaChannel.lon(MOORING.KwakshuaChannel.lon<0) + 360;

%% Dabob
% MOORING.Dabob.name = name;
% MOORING.Dabob.lat  = lat;
% MOORING.Dabob.lon  = lon;
% MOORING.Dabob.date = date;
% MOORING.Dabob.time = time;
% MOORING.Dabob.sst  = sst;
% MOORING.Dabob.sal  = sal;
% MOORING.Dabob.fCO2SW  = fCO2SW;
% MOORING.Dabob.fCO2Air = fCO2Air;
% idx = find(strcmp('Dabob',structnames));
% for n = 1:size(idx,1)
%     HL = 0; if idx(n)==24 || idx(n)==28; HL = 4; end
%     MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
%         readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
%     idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
%     MOORING.Dabob.name    = [MOORING.Dabob.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
%     MOORING.Dabob.lat     = [MOORING.Dabob.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
%     MOORING.Dabob.lon     = [MOORING.Dabob.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
%     MOORING.Dabob.date    = [MOORING.Dabob.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
%     %MOORING.Dabob.time   = [MOORING.Dabob.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
%     MOORING.Dabob.sst     = [MOORING.Dabob.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
%     MOORING.Dabob.sal     = [MOORING.Dabob.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
%     MOORING.Dabob.fCO2SW  = [MOORING.Dabob.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
%     MOORING.Dabob.fCO2Air = [MOORING.Dabob.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
% end
% MOORING.Dabob.date(MOORING.Dabob.date(:,1)<1000)=MOORING.Dabob.date(MOORING.Dabob.date(:,1)<1000)+2000;
% MOORING.Dabob.sst(MOORING.Dabob.sst==-999) = NaN;
% MOORING.Dabob.sal(MOORING.Dabob.sal==-999) = NaN;
% MOORING.Dabob.fCO2SW(MOORING.Dabob.fCO2SW==-999) = NaN;
% MOORING.Dabob.fCO2Air(MOORING.Dabob.fCO2Air==-999) = NaN;
% MOORING.Dabob.lon(MOORING.Dabob.lon<0) = MOORING.Dabob.lon(MOORING.Dabob.lon<0) + 360;

%% Twanoh
% MOORING.Twanoh.name = name;
% MOORING.Twanoh.lat  = lat;
% MOORING.Twanoh.lon  = lon;
% MOORING.Twanoh.date = date;
% MOORING.Twanoh.time = time;
% MOORING.Twanoh.sst  = sst;
% MOORING.Twanoh.sal  = sal;
% MOORING.Twanoh.fCO2SW  = fCO2SW;
% MOORING.Twanoh.fCO2Air = fCO2Air;
% idx = find(strcmp('Twanoh',structnames));
% for n = 1:size(idx,1)
%     HL = 0; if idx(n)==68 || idx(n)==69 || idx(n)==70; HL = 4; end
%     MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
%         readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
%     if n == 5 || n==6 || n==7
%         idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
%     else
%         idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
%     end
%     MOORING.Twanoh.name    = [MOORING.Twanoh.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
%     MOORING.Twanoh.lat     = [MOORING.Twanoh.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
%     MOORING.Twanoh.lon     = [MOORING.Twanoh.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
%     MOORING.Twanoh.date    = [MOORING.Twanoh.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
%     %MOORING.Twanoh.time   = [MOORING.Twanoh.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
%     MOORING.Twanoh.sst     = [MOORING.Twanoh.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
%     MOORING.Twanoh.sal     = [MOORING.Twanoh.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
%     MOORING.Twanoh.fCO2SW  = [MOORING.Twanoh.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
%     MOORING.Twanoh.fCO2Air = [MOORING.Twanoh.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
% end
% MOORING.Twanoh.date(MOORING.Twanoh.date(:,1)<1000)=MOORING.Twanoh.date(MOORING.Twanoh.date(:,1)<1000)+2000;
% MOORING.Twanoh.sst(MOORING.Twanoh.sst==-999) = NaN;
% MOORING.Twanoh.sal(MOORING.Twanoh.sal==-999) = NaN;
% MOORING.Twanoh.fCO2SW(MOORING.Twanoh.fCO2SW==-999) = NaN;
% MOORING.Twanoh.fCO2Air(MOORING.Twanoh.fCO2Air==-999) = NaN;
% MOORING.Twanoh.lon(MOORING.Twanoh.lon<0) = MOORING.Twanoh.lon(MOORING.Twanoh.lon<0) + 360;

%% Papa
% MOORING.Papa.name = name;
% MOORING.Papa.lat  = lat;
% MOORING.Papa.lon  = lon;
% MOORING.Papa.date = date;
% MOORING.Papa.time = time;
% MOORING.Papa.sst  = sst;
% MOORING.Papa.sal  = sal;
% MOORING.Papa.fCO2SW  = fCO2SW;
% MOORING.Papa.fCO2Air = fCO2Air;
% idx = find(strcmp('Papa',structnames));
% for n = 1:size(idx,1)
%     HL = 4; if n==5 || n==6; HL = 0; end
%     MOORING.(extractBefore(fnames2{idx(n)},'.')) = ...
%         readtable(strcat(filePath,fnames{idx(n)}),'HeaderLines',HL,'VariableNamingRule','modify');
%     if n == 5 || n==6
%     idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).QF == 2;
%     else
%     idxQC = MOORING.(extractBefore(fnames2{idx(n)},'.')).CO2SWQF == 2;
%     end
%     MOORING.Papa.name    = [MOORING.Papa.name;MOORING.(extractBefore(fnames2{idx(n)},'.')).MooringName(idxQC)];
%     MOORING.Papa.lat     = [MOORING.Papa.lat;MOORING.(extractBefore(fnames2{idx(n)},'.')).Latitude(idxQC)];
%     MOORING.Papa.lon     = [MOORING.Papa.lon;MOORING.(extractBefore(fnames2{idx(n)},'.')).Longitude(idxQC)];
%     MOORING.Papa.date    = [MOORING.Papa.date;datevec(MOORING.(extractBefore(fnames2{idx(n)},'.')).Date(idxQC))];
%     %MOORING.Papa.time   = [MOORING.Papa.time;MOORING.(extractBefore(fnames2{idx(n)},'.')).Time(idxQC)];
%     MOORING.Papa.sst     = [MOORING.Papa.sst;MOORING.(extractBefore(fnames2{idx(n)},'.')).SST_C_(idxQC)];
%     MOORING.Papa.sal     = [MOORING.Papa.sal;MOORING.(extractBefore(fnames2{idx(n)},'.')).Salinity(idxQC)];
%     MOORING.Papa.fCO2SW  = [MOORING.Papa.fCO2SW;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2SW_sat_Uatm(idxQC)];
%     MOORING.Papa.fCO2Air = [MOORING.Papa.fCO2Air;MOORING.(extractBefore(fnames2{idx(n)},'.')).fCO2Air_sat_Uatm(idxQC)];
% end
% MOORING.Papa.date(MOORING.Papa.date(:,1)<1000)=MOORING.Papa.date(MOORING.Papa.date(:,1)<1000)+2000;
% MOORING.Papa.sst(MOORING.Papa.sst==-999) = NaN;
% MOORING.Papa.sal(MOORING.Papa.sal==-999) = NaN;
% MOORING.Papa.fCO2SW(MOORING.Papa.fCO2SW==-999) = NaN;
% MOORING.Papa.fCO2Air(MOORING.Papa.fCO2Air==-999) = NaN;
% MOORING.Papa.lon(MOORING.Papa.lon<0) = MOORING.Papa.lon(MOORING.Papa.lon<0) + 360;

%% Remove extraneous fields from MOORING structure
fields = fieldnames(MOORING);
idx = contains(fields,'_');
idx = idx & ~strcmp(fields,'CB_06');
MOORING = rmfield(MOORING,fields(idx));

%% Clean up Workspace
clear date fileinfo fnames fnames2 HL idx lat lon n name filePath
clear fields fCO2Air fCO2SW sal sst structnames time idxQC hr mn

%% Define mooring names
moornames = fieldnames(MOORING);

%% Bin mooring fCO2s into monthly values
for n = 1:numel(moornames)
    years = unique(MOORING.(moornames{n}).date(:,1));
    years = sort(years);
    MOORING.(moornames{n}).sst_monthly = nan(length(years)*12,2);
    MOORING.(moornames{n}).sal_monthly = nan(length(years)*12,2);
    MOORING.(moornames{n}).fCO2SW_monthly = nan(length(years)*12,2);
    MOORING.(moornames{n}).fCO2Air_monthly = nan(length(years)*12,2);
    for y = 1:length(years)
        for m = 1:12
            idx = MOORING.(moornames{n}).date(:,1)==years(y) & ...
                MOORING.(moornames{n}).date(:,2)==m;
            MOORING.(moornames{n}).datetime_monthly((y-1)*12+m,1) = ...
                datenum(years(y),m,15);
            MOORING.(moornames{n}).sst_monthly((y-1)*12+m,1) = ...
                mean(MOORING.(moornames{n}).sst(idx),'omitnan');
            MOORING.(moornames{n}).sal_monthly((y-1)*12+m,1) = ...
                mean(MOORING.(moornames{n}).sal(idx),'omitnan');
            MOORING.(moornames{n}).fCO2SW_monthly((y-1)*12+m,1) = ...
                mean(MOORING.(moornames{n}).fCO2SW(idx),'omitnan');
            MOORING.(moornames{n}).fCO2Air_monthly((y-1)*12+m,1) = ...
                mean(MOORING.(moornames{n}).fCO2Air(idx),'omitnan');
            MOORING.(moornames{n}).sst_monthly((y-1)*12+m,2) = ...
                std(MOORING.(moornames{n}).sst(idx),[],'omitnan');
            MOORING.(moornames{n}).sal_monthly((y-1)*12+m,2) = ...
                std(MOORING.(moornames{n}).sal(idx),[],'omitnan');
            MOORING.(moornames{n}).fCO2SW_monthly((y-1)*12+m,2) = ...
                std(MOORING.(moornames{n}).fCO2SW(idx),[],'omitnan');
            MOORING.(moornames{n}).fCO2Air_monthly((y-1)*12+m,2) = ...
                std(MOORING.(moornames{n}).fCO2Air(idx),[],'omitnan');
        end
    end
end
clear m n

%% Bin mooring fCO2s into climatological values
for n = 1:numel(moornames)
    MOORING.(moornames{n}).sst_clim = nan(12,2);
    MOORING.(moornames{n}).sal_clim = nan(12,2);
    MOORING.(moornames{n}).fCO2SW_clim = nan(12,2);
    MOORING.(moornames{n}).fCO2Air_clim = nan(12,2);
    % calculate climatological values
    for m = 1:12
        idx = MOORING.(moornames{n}).date(:,2)==m;
        MOORING.(moornames{n}).sst_clim(m,1) = ...
            mean(MOORING.(moornames{n}).sst(idx),'omitnan');
        MOORING.(moornames{n}).sal_clim(m,1) = ...
            mean(MOORING.(moornames{n}).sal(idx),'omitnan');
        MOORING.(moornames{n}).fCO2SW_clim(m,1) = ...
            mean(MOORING.(moornames{n}).fCO2SW(idx),'omitnan');
        MOORING.(moornames{n}).fCO2Air_clim(m,1) = ...
            mean(MOORING.(moornames{n}).fCO2Air(idx),'omitnan');
        MOORING.(moornames{n}).sst_clim(m,2) = ...
            std(MOORING.(moornames{n}).sst(idx),[],'omitnan');
        MOORING.(moornames{n}).sal_clim(m,2) = ...
            std(MOORING.(moornames{n}).sal(idx),[],'omitnan');
        MOORING.(moornames{n}).fCO2SW_clim(m,2) = ...
            std(MOORING.(moornames{n}).fCO2SW(idx),[],'omitnan');
        MOORING.(moornames{n}).fCO2Air_clim(m,2) = ...
            std(MOORING.(moornames{n}).fCO2Air(idx),[],'omitnan');
    end
end
clear m n

%% Calculate amplitudes, trends, and annual means
for n = 1:numel(moornames)
    % index to where fCO2 is available
    idx = ~isnan(MOORING.(moornames{n}).fCO2SW_monthly(:,1));
    t = MOORING.(moornames{n}).datetime_monthly;
    t_idx = MOORING.(moornames{n}).datetime_monthly(idx);
    dt_idx = t_idx-t(1);
    fCO2 = MOORING.(moornames{n}).fCO2SW_monthly(idx,1);
    % fit regression model
    [~,~,x] = leastsq2(dt_idx,fCO2,0,2,[365.25 365.25/2]);
    % calculate fCO2 from model
    dt = t-t(1);
    freqp1 = 2*pi/365.25;
    freqp2 = 2*pi/(365.25/2);
    yf = x(1) + x(2)*dt + x(3)*cos(freqp1*dt) + x(4)*sin(freqp1*dt) + ...
        x(5)*cos(freqp2*dt) + x(6)*sin(freqp2*dt);
    % plot
%     figure; hold on
%     scatter(t_idx,fCO2);
%     plot(t,yf);
    % compute annual mean
    MOORING.(moornames{n}).fCO2SW_annmean = mean(yf);
    % compute amplitude
    yf_no_tr = x(1) + x(3)*cos(freqp1*dt) + x(4)*sin(freqp1*dt) + ...
        x(5)*cos(freqp2*dt) + x(6)*sin(freqp2*dt);
    MOORING.(moornames{n}).fCO2SW_amplitude = max(yf_no_tr) - min(yf_no_tr);
    % compute trend
    if length(years)>6
        MOORING.(moornames{n}).fCO2SW_trend = x(2)*365.25;
    else
        MOORING.(moornames{n}).fCO2SW_trend = NaN;
    end
end