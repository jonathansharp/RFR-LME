% Load SOCAT data

function load_socat(vrs)

% display status
disp('Loading SOCAT data');

% describe data format and parts to skip
fmt_data = ['%s %*s %*s %s %s %s %s %*s %*s %*s %s %s %*s %s %s %*s %*s %*s ' ...
    '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s %*s %s \n'];

% process header lines to find where data starts
ff = fopen([vrs '/' vrs '.tsv'],'r'); % open .tsv file
header_lines = textscan(ff,'%s',10000,'delimiter','\n'); % scan first 10,000 lines
hdr_in1 = strfind(header_lines{1}, 'Expocode'); hdr_in1=~cellfun('isempty', hdr_in1);
hdr_in2 = strfind(header_lines{1}, 'SST'); hdr_in2=~cellfun('isempty', hdr_in2);
hdrline=intersect(find(hdr_in1), find(hdr_in2)); % number of header lines before data starts
clear header_lines
fclose(ff);

% import data
ff = fopen([vrs '/' vrs '.tsv'],'r'); % open .tsv file
junk = textscan(ff,'%s',hdrline-1,'delimiter','\n'); clear junk;
headers = textscan(ff,fmt_data,1,'delimiter','\t');
headers = cellfun(@eraseBetween,headers,repmat({' ['},1,length(headers)),...
    repmat({']'},1,length(headers)),repmat({'Boundaries'},1,length(headers)),...
    repmat({'inclusive'},1,length(headers)));
for h = 1:length(headers); data.(headers{h}) = []; end % pre-allocate
lat_idx = find(cellfun(@contains,headers,repmat({'latitude'},1,11)));
lon_idx = find(cellfun(@contains,headers,repmat({'longitude'},1,11)));
n=5000000; c = 1;
while n == 5000000 % when we get to the end of the file, the last scan will have fewer than 5,000,000 values
    data_scn = textscan(ff,fmt_data,5000000,'delimiter','\t'); % scan 5,000,000 lines
    lat = str2double(data_scn{1,lat_idx}); % extract latitude
    lon = str2double(data_scn{1,lon_idx}); % extract longitude
    idx = lat > -18 & lat < 82 & lon > 140 & lon < 302; % range of interest
    if sum(idx) > 0 % if any obs fall within the spatial range
        for h = 1:length(headers)
            % pull data from scanned lines
            if h <= 2
                data.(headers{h}) = [data.(headers{h});data_scn{1,h}(idx)];
            else
                data.(headers{h}) = [data.(headers{h});str2double(data_scn{1,h}(idx))];
            end
        end
    end
    if mod(c,1) == 0; disp([num2str(c*n) ' lines scanned']); end
    n = length(data_scn{1}); c = c + 1;
    clear data_scn
end
disp('all lines scanned');
fclose(ff);


% display status
disp('Processing SOCAT data');

% remove observations before 1998
idxyr = data.yr >= 1998;
vars = fieldnames(data);
for v = 1:numel(vars)
    tempvar = data.(string(vars(v)));
    data.(string(vars(v))) = tempvar(idxyr);
end

% determine months since 1 Jan 1998
data.month_since_1998 = (data.yr-1998).*12 + data.mon;

% remove observations with flags other than 2 and A/B/C/D
idxflag = data.fCO2rec_flag == 2 & ...
    (strcmp(data.QC_Flag,'A') | strcmp(data.QC_Flag,'B') | ...
     strcmp(data.QC_Flag,'C') | strcmp(data.QC_Flag,'D'));
vars = fieldnames(data);
for v = 1:numel(vars)
    tempvar = data.(char(vars(v)));
    data.(char(vars(v))) = tempvar(idxflag);
end

% obtain unique integers for each expocode
data.cruise = nan(size(data.Expocode));
cruiselist = unique(data.Expocode);
for c = 1:numel(cruiselist)
    idx = strcmp(cruiselist(c),data.Expocode);
    data.cruise(idx) = c;
end

% visualize number of observations
time = datetime(data.yr,data.mon,data.day);
figure('visible','on');
histogram(time);
ylabel('Number of observations');
xlabel('Year');
clear time

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/fCO2_hist_NA.png');
close

% save SOCAT structure
if ~isfolder('Data'); mkdir('Data'); end
save(['Data/' vrs '_structure'],'data','-v7.3');

end