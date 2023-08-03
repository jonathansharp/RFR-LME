% Run scripts to evaluate US LME OA indicators

%% load CODAP data
load('CODAP_NA/CODAP_NA_v2020_G2format.mat')
idx = C1.pressure <= 10 & C1.tco2f == 2 & C1.talkf == 2 & ...
    C1.salinityf == 2 & C1.silicatef == 2 & C1.phosphatef == 2;
CODAP.lat = C1.latitude(idx);
CODAP.lon = C1.longitude(idx);
CODAP.lon = convert_lon(CODAP.lon);
CODAP.time = datenum(C1.year(idx),C1.month(idx),C1.day(idx));
CODAP.prs = C1.pressure(idx);
CODAP.sal = C1.salinity(idx);
CODAP.tmp = C1.temperature(idx);
CODAP.sil = C1.silicate(idx);
CODAP.phos = C1.phosphate(idx);
CODAP.DIC = C1.tco2(idx);
CODAP.TA = C1.talk(idx);
CODAP.pH = C1.phtsinsitutp(idx);
CODAP.OmA = C1.aragonite(idx);
CODAP.OmC = C1.calcite(idx);
CODAP.RF = C1.revelle(idx);
clear C1 idx Headers Units

%% load GLODAP data
load('GLODAPv2.2022/GLODAPv2.2022_Merged_Master_File.mat')
idx = G2pressure <= 10 & G2tco2f == 2 & G2talkf == 2 & ...
    G2salinityf == 2 & G2silicatef == 2 & G2phosphatef == 2;
GLODAP.lat = G2latitude(idx);
GLODAP.lon = G2longitude(idx);
GLODAP.lon = convert_lon(GLODAP.lon);
GLODAP.time = datenum(G2year(idx),G2month(idx),G2day(idx));
GLODAP.prs = G2pressure(idx);
GLODAP.sal = G2salinity(idx);
GLODAP.tmp = G2temperature(idx);
GLODAP.sil = G2silicate(idx);
GLODAP.phos = G2phosphate(idx);
GLODAP.DIC = G2tco2(idx);
GLODAP.TA = G2talk(idx);
GLODAP.pH = G2phtsinsitutp(idx);
clear idx G2* expocode expocodeno

%% remove CODAP data outside LMEs
% determine LME index
define_regions_eiwg
idx_tmp = nan(length(CODAP.lat),length(region));
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_tmp_CODAP(:,n) = inpolygon(CODAP.lon,CODAP.lat,tmp_lon,tmp_lat);
    idx_tmp_GLODAP(:,n) = inpolygon(GLODAP.lon,GLODAP.lat,tmp_lon,tmp_lat);
end
idx_CODAP = any(idx_tmp_CODAP,2);
idx_GLODAP = any(idx_tmp_GLODAP,2);
% remove outside data (CODAP)
vars = fieldnames(CODAP);
for v = 1:length(vars)
    CODAP.(vars{v}) = CODAP.(vars{v})(idx_CODAP);
end
% remove outside data (GLODAP)
vars = fieldnames(GLODAP);
for v = 1:length(vars)
    GLODAP.(vars{v}) = GLODAP.(vars{v})(idx_GLODAP);
end
clear idx_CODAP idx_tmp_CODAP idx_GLODAP idx_tmp_GLODAP vars v tmp_lon tmp_lat

%% calculate fCO2, H, and CO3 (CODAP)
carb = CO2SYS(CODAP.TA,CODAP.DIC,1,2,CODAP.sal,CODAP.tmp,NaN,CODAP.prs,...
    NaN,CODAP.sil,CODAP.phos,0,0,1,10,1,2,2);
CODAP.pco2 = carb(:,4);
CODAP.fco2 = carb(:,5);
CODAP.H = carb(:,15).*10^3;
CODAP.CO3 = carb(:,7);
clear carb

%% calculate fCO2, H, CO3, Omega, and RF (GLODAP)
carb = CO2SYS(GLODAP.TA,GLODAP.DIC,1,2,GLODAP.sal,GLODAP.tmp,NaN,GLODAP.prs,...
    NaN,GLODAP.sil,GLODAP.phos,0,0,1,10,1,2,2);
GLODAP.pco2 = carb(:,5);
GLODAP.fco2 = carb(:,5);
GLODAP.H = carb(:,15).*10^3;
GLODAP.CO3 = carb(:,7);
GLODAP.OmA = carb(:,18);
GLODAP.OmC = carb(:,17);
GLODAP.RF = carb(:,16);
clear carb

%% load US LME data
date = '02-Aug-2023';
LME_RFR.lat = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Lat');
LME_RFR.lon = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Lon');
LME_RFR.time = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'Time');
LME_RFR.pco2 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'pCO2');
LME_RFR.fco2 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'fCO2');
LME_RFR.TA = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'TA');
LME_RFR.DIC = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'DIC');
LME_RFR.pH = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'pH');
LME_RFR.OmA = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'OmA');
LME_RFR.OmC = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'OmC');
LME_RFR.H = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'H');
LME_RFR.CO3 = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'CO3');
LME_RFR.RF = ncread(['Data/US_LME_RFR_Inds_' date '.nc'],'RF');

%% plot data locations
define_regions_eiwg
% initialize figure
figure('visible','on'); box on; hold on;
worldmap([-18 82],[140 302]);
setm(gca,'MapProjection','robinson','MLabelParallel','south');
set(gcf,'position',[100 100 900 600]);
set(gca,'fontsize',16);
% figure properties
c=colorbar('location','southoutside');
colormap(parula(18));
caxis([295 475]);
c.TickLength = 0;
c.Label.String = 'Sea Surface {\itp}CO_{2}';
cbarrow;
% plot regions
for n = 1:length(region)
     load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    z = mean(OAI_grid.(region{n}).pCO2,3,'omitnan')';
    contourfm(OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        z,295:10:475,'LineStyle','none');
    clear vars_grid z
end
% plot borders around regions
for n = 1:length(region)
    if n <= 11
        tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    else
        tmp_lon = lme_shape(lme_idx.(region{n})).X';
    end
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    plotm(tmp_lat,tmp_lon,'k','linewidth',1);
end
% plot land
plot_land('map');
mlabel off
% add CODAP to plot
scatterm(CODAP.lat,CODAP.lon,5,'ro')
% add GLODAP to plot
scatterm(GLODAP.lat,GLODAP.lon,5,'g.')
% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,'Figures/CODAP_GLODAP_eval_map.png');

%% variable information
edges = {1900:10:2300;300:5:500;300:5:500;2000:10:2400;7.8:0.0125:8.3;...
    1:0.1:5;1:0.1:5;7:0.2:14;100:2.5:250;9:0.1:17};
var_type = {'DIC' 'fco2' 'pco2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}_{CO2}' '{\itp}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 1 3 2 2 1 1 2];

%% co-locate each CODAP point with an LME-RFR grid cell
% pre-allocate
for v = 1:length(var_type)
    CODAP.([var_type{v} '_grid']) = nan(size(CODAP.(var_type{v})));
end
% identify and fill gridded values
for n = 1:length(CODAP.lat)
    idx_lon = find(abs(LME_RFR.lon - CODAP.lon(n)) == ...
        min(abs(LME_RFR.lon - CODAP.lon(n))));
    idx_lat = find(abs(LME_RFR.lat - CODAP.lat(n)) == ...
        min(abs(LME_RFR.lat - CODAP.lat(n))));
    idx_time = find(abs(LME_RFR.time - CODAP.time(n)) == ...
        min(abs(LME_RFR.time - CODAP.time(n))));
    for v = 2:length(var_type)
        CODAP.([var_type{v} '_grid'])(n) = ...
            LME_RFR.(var_type{v})(idx_lon(1),idx_lat(1),idx_time(1));
    end
end
% calculate differences
for v = 2:length(var_type)
    CODAP.([var_type{v} '_del']) = ...
        CODAP.(var_type{v}) - CODAP.([var_type{v} '_grid']);
end

%% co-locate each GLODAP point with an LME-RFR grid cell
% pre-allocate
for v = 1:length(var_type)
    GLODAP.([var_type{v} '_grid']) = nan(size(GLODAP.(var_type{v})));
end
% identify and fill gridded values
for n = 1:length(GLODAP.lat)
    idx_lon = find(abs(LME_RFR.lon - GLODAP.lon(n)) == ...
        min(abs(LME_RFR.lon - GLODAP.lon(n))));
    idx_lat = find(abs(LME_RFR.lat - GLODAP.lat(n)) == ...
        min(abs(LME_RFR.lat - GLODAP.lat(n))));
    idx_time = find(abs(LME_RFR.time - GLODAP.time(n)) == ...
        min(abs(LME_RFR.time - GLODAP.time(n))));
    for v = 2:length(var_type)
        GLODAP.([var_type{v} '_grid'])(n) = ...
            LME_RFR.(var_type{v})(idx_lon(1),idx_lat(1),idx_time(1));
    end
end
% eliminate NaN and -999
for v = 1:length(var_type)
    idx = ~isnan(GLODAP.(var_type{v})) & GLODAP.(var_type{v}) > 0;
    GLODAP.(var_type{v}) = GLODAP.(var_type{v})(idx);
    GLODAP.([var_type{v} '_grid']) = GLODAP.([var_type{v} '_grid'])(idx);
end
% calculate differences
for v = 1:length(var_type)
    GLODAP.([var_type{v} '_del']) = ...
        GLODAP.(var_type{v}) - GLODAP.([var_type{v} '_grid']);
end

%% plot figures
% for v = 2:length(var_type)
%     plot_delta_eval(edges{v},CODAP.(var_type{v}),...
%         CODAP.([var_type{v} '_grid']),var_type{v},var_lab{v},units{v},rounder(v),'CODAP')
% end

%% TA
% CODAP
figure; hold on;
title('TA');
scatter(CODAP.TA,CODAP.TA_grid,'k.');
plot([0 2500],[0 2500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([1500 2500]);
ylim([1500 2500]);
text(1500,2400,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.TA_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.TA_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_TA.png');
% GLODAP
figure; hold on;
title('TA');
scatter(GLODAP.TA,GLODAP.TA_grid,'k.');
plot([1500 2500],[1500 2500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([1500 2500]);
ylim([1500 2500]);
text(1500,2400,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.TA_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.TA_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_TA.png');

%% fCO2
% CODAP
figure; hold on;
title('fCO2');
scatter(CODAP.fco2,CODAP.fco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_fCO2.png');
close
% GLODAP
figure; hold on;
title('fCO2');
scatter(GLODAP.fco2,GLODAP.fco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_fCO2.png');
close
% Both
figure; hold on;
scatter(GLODAP.fco2,GLODAP.fco2_grid,'g.');
scatter(CODAP.fco2,CODAP.fco2_grid,'r.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP/CODAP {\itf}CO_{2}');
ylabel('LME Grid {\itf}CO_{2}');
xlim([0 1200]);
ylim([0 1200]);
text(200,1100,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.fco2_del,[],'omitnan'),1))]);
text(200,900,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.fco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.fco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_fCO2.png');
close
% Both histogram
figure; hold on;
histogram(GLODAP.fco2_del);
histogram(CODAP.fco2_del);
xlabel('{\itf}_{CO2(discrete)} - {\itf}_{CO2(RFR)}');
ylabel('Counts');
xlim([-500 1000]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_fCO2.png');
close

%% pCO2
% CODAP
figure; hold on;
title('pCO2');
scatter(CODAP.pco2,CODAP.pco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/CODAP_eval_pCO2.png');
close
% GLODAP
figure; hold on;
title('pCO2');
scatter(GLODAP.pco2,GLODAP.pco2_grid,'k.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([0 1200]);
ylim([0 1200]);
text(200,1000,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_eval_pCO2.png');
close
% Both
figure; hold on;
scatter(GLODAP.pco2,GLODAP.pco2_grid,'g.');
scatter(CODAP.pco2,CODAP.pco2_grid,'r.');
plot([0 1500],[0 1500],'k-');
xlabel('GLODAP/CODAP {\itp}CO_{2}');
ylabel('LME Grid {\itp}CO_{2}');
xlim([0 1200]);
ylim([0 1200]);
text(200,1100,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pco2_del,[],'omitnan'),1))]);
text(200,900,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pco2_del,'omitnan'),1)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pco2_del,[],'omitnan'),1))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_pCO2.png');
close
% Both histogram
figure; hold on;
histogram(GLODAP.pco2_del);
histogram(CODAP.pco2_del);
plot([0 0],[0 800],'--k','linewidth',2);
xlabel('{\itp}CO_{2(discrete)} - {\itp}CO_{2(RFR)}');
ylabel('Counts');
xlim([-500 1000]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_pCO2.png');
close

%% pH
% CODAP
figure; hold on;
title('pH');
scatter(CODAP.pH,CODAP.pH_grid,'k.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlabel('CODAP');
ylabel('LME Grid');
xlim([7.4 8.5]);
ylim([7.4 8.5]);
text(7.5,8.4,['CODAP ' char(45) ' RFR-LME = ' num2str(round(mean(CODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/CODAP_eval_pH.png');
close
% GLODAP
figure; hold on;
title('pH');
scatter(GLODAP.pH,GLODAP.pH_grid,'k.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
xlim([7.4 8.6]);
ylim([7.4 8.6]);
text(7.5,8.4,['GLODAP ' char(45) ' RFR-LME = ' num2str(round(mean(GLODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_eval_pH.png');
close
% Both scatter
figure; hold on;
scatter(GLODAP.pH,GLODAP.pH_grid,'g.');
scatter(CODAP.pH,CODAP.pH_grid,'r.');
plot([7.3 8.7],[7.3 8.7],'k-');
xlim([7.4 8.6]);
ylim([7.4 8.6]);
xlabel('GLODAP/CODAP pH_{T}');
ylabel('LME Grid pH_{T}');
text(7.5,8.5,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.pH_del,[],'omitnan'),3))]);
text(7.5,8.42,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.pH_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.pH_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_pH.png');
% Both histogram
figure; hold on;
histogram(GLODAP.pH_del);
histogram(CODAP.pH_del);
plot([0 0],[0 900],'--k','linewidth',2);
xlabel('pH_{T(discrete)} - pH_{T(RFR)}');
ylabel('Counts');
xlim([-1 1]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_pH.png');
close

%% OmA
% CODAP
figure; hold on;
title('\Omega_{A}');
scatter(CODAP.OmA,CODAP.OmA_grid,'k.');
plot([0 5],[0 5],'k-');
xlabel('CODAP');
ylabel('LME Grid');
text(0.5,4.5,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/CODAP_eval_OmA.png');
close
% GLODAP
figure; hold on;
title('\Omega_{A}');
scatter(GLODAP.OmA,GLODAP.OmA_grid,'k.');
plot([0 5],[0 5],'k-');
xlabel('GLODAP');
ylabel('LME Grid');
text(0.5,4.5,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_eval_OmA.png');
close
% Both scatter
figure; hold on;
scatter(GLODAP.OmA,GLODAP.OmA_grid,'g.');
scatter(CODAP.OmA,CODAP.OmA_grid,'r.');
plot([0 5],[0 5],'k-');
ylim([0 5]); xlim([0 5]);
xlabel('GLODAP/CODAP \Omega_{A}');
ylabel('LME Grid \Omega_{A}');
text(0.5,4.7,['GLODAP ' char(45) ' LME = ' num2str(round(mean(GLODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(GLODAP.OmA_del,[],'omitnan'),3))]);
text(0.5,4.4,['CODAP ' char(45) ' LME = ' num2str(round(mean(CODAP.OmA_del,'omitnan'),3)),...
    ' ' char(177) ' ' num2str(round(std(CODAP.OmA_del,[],'omitnan'),3))]);
exportgraphics(gcf,'Figures/GLODAP_CODAP_eval_OmA.png');
% Both histogram
figure; hold on;
histogram(GLODAP.OmA_del);
histogram(CODAP.OmA_del);
plot([0 0],[0 450],'--k','linewidth',2);
xlabel('\Omega_{A(discrete)} - \Omega_{A(RFR)}');
ylabel('Counts');
xlim([-2.5 2.5]);
legend({'GLODAP' 'CODAP'});
exportgraphics(gcf,'Figures/hist_OmA.png');
close
