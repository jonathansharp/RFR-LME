%% compare US-LME-RFR to time series data

% load CalCOFI data
calcofi = csvread('raid/Data/CalCOFI/194903-201911_Bottle.csv',[R1 C1 R2 C2]);

% load GLODAP data
load('raid/Data/GLODAPv2.2022/GLODAPv2.2022_Merged_Master_File.mat');

