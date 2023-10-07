%% plot station locations for SOCAT, GLODAP, CODAP

% load GLODAP data
load('GLODAPv2.2022/GLODAPv2.2022_Merged_Master_File.mat')
GLODAP.lat = G2latitude;
GLODAP.lon = G2longitude;
GLODAP.idx = G2talkf == 2 & G2tco2f == 2;
clearvars -except GLODAP
% load CODAP data
load('CODAP_NA/CODAP_NA_v2020_G2format.mat')
CODAP.lat = C1.latitude;
CODAP.lon = C1.longitude;
CODAP.idx = C1.talkf == 2 & C1.tco2f == 2;
clearvars -except GLODAP CODAP

% define regions
define_regions
% plot station locations for each LME
for n = 1:length(region)
    % load regional SOCAT
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    % initialize figure
    figure;
    worldmap([min(Vars_array.(region{n}).X_mod(:,2)) ...
              max(Vars_array.(region{n}).X_mod(:,2))],...
             [min(Vars_array.(region{n}).X_mod(:,1)) ...
              max(Vars_array.(region{n}).X_mod(:,1))]);
    % plot SOCAT data locations
    scatterm(Vars_array.(region{n}).X_mod(:,2),...
        Vars_array.(region{n}).X_mod(:,1),'k.');
    % plot GLODAP data locations
    scatterm(GLODAP.lat(GLODAP.idx),GLODAP.lon(GLODAP.idx),'bo');
    % plot CODAP data locations
    scatterm(CODAP.lat(CODAP.idx),CODAP.lon(CODAP.idx),'r.');
    % plot land
    plot_land('map');
    % plot borders around regions
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
        plotm(tmp_lat,tmp_lon,'k','linewidth',1);
    % save figure
    exportgraphics(gcf,['stn_locs/' region{n} '.png']);
    close
    % clean up
    clear Vars_array
end
