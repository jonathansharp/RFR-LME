try

% Script to create LME-RFR indicators
vrs = 'SOCATv2024';
dpath = '/raid/Data/';

% load_socat(vrs);
% grid_socat(vrs);
load_vars(vrs,dpath);
% extract_lme(vrs);

catch
    exit
end