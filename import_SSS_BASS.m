%% Import BASS salinity
file = fopen([path '/Data/BASS_V0.Z_MON_1DEG.lnx.B201001.txt']);
ECCO_SSS.SSS = fread(file);