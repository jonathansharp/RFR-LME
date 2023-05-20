% Plot individual OA indicator time series
% 
% Written by J.D. Sharp: 2/20/23
% Last updated by J.D. Sharp: 2/20/23
% 

% variable type
var_num = 1;
% region number
n = 9;

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' ...
           'HI' 'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};
% variable information
var_type = {'DIC' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 3 2 2 1 1 2];
% initialize figure
figure('Visible','on'); hold on;
set(gcf,'position',[10 10 900 400]);
% load estimated OA grid
load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
% calculate area-weighted time series
OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
for t = 1:OAI_grid.(region{n}).dim.z
    OAI_grid.(region{n}).var_dom_mean(t) = ...
        squeeze(sum(sum(OAI_grid.(region{n}).(var_type{var_num})(:,:,t).*...
            area_weights,1,'omitnan'),2,'omitnan'))./...
            squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
end
% scale H to nanomoles
if strcmp(var_type{var_num},'H')
    OAI_grid.(region{n}).var_dom_mean = (10^9).*OAI_grid.(region{n}).var_dom_mean;
end
% re-calculate time
time = datenum([OAI_grid.(region{n}).year ...
                OAI_grid.(region{n}).month_of_year ...
                repmat(15,OAI_grid.(region{n}).dim.z,1)]);
% calculate annual means
for y = 1:length(OAI_grid.(region{n}).year)/12
    year(y) = mean(time((y-1)*12+1:(y-1)*12+12));
    var_ann(y) = mean(OAI_grid.(region{n}).var_dom_mean((y-1)*12+1:(y-1)*12+12));
end
% calculate trend
[yf,yr,x] = ...
    leastsq2(OAI_grid.(region{n}).month,...
    OAI_grid.(region{n}).var_dom_mean,0,2,[6 12]);

% calculate trend over most recent 5 years
[yf5,yr5,x5,err5] = ...
    leastsq2(OAI_grid.(region{n}).month(end-12*5:end),...
    OAI_grid.(region{n}).var_dom_mean(end-12*5:end),0,2,[6 12]);
pos_signif_5 = x5(2)-err5(2) > 0;
neg_signif_5 = x5(2)+err5(2) < 0;
% calculate mean over most recent 5 years
meantot = mean(OAI_grid.(region{n}).var_dom_mean);
conf90 = std(OAI_grid.(region{n}).var_dom_mean)*1.645;
mean5 = mean(OAI_grid.(region{n}).var_dom_mean(end-12*5:end));
pos_mean = mean5 > meantot+conf90;
neg_mean = mean5 < meantot-conf90;

% determine average modelled climatology
clim = nan(12,1);
for m = 1:12
    clim(m) = mean(yf(m:12:end));
end
% calculate amplitude
amp = max(clim) - min(clim);

% plot fCO2 time series
plot(time,OAI_grid.(region{n}).var_dom_mean,'k-','linewidth',1);
plot(year,var_ann,'r-','linewidth',3);
scatter(year,var_ann,200,'k.');
datetick('x');
xlim([datenum([1998 1 1]) datenum([2022 1 1])]);   
% plot means
plot([time(1) time(end-12*5)],[meantot meantot],'k--','linewidth',2);
plot([time(end-12*5) time(end)],[mean5 mean5],'k--','linewidth',2);
% add text to plot
title([reg_lab{n} ' | \mu = ' num2str(round(mean(OAI_grid.(region{n}).var_dom_mean),rounder(var_num))) ...
    ' ' units{var_num} ' | ' num2str(round(x(2)*12,rounder(var_num))) ...
    ' ' units{var_num} ' yr^{-1} | Amp. = ' ...
    num2str(round(amp,rounder(var_num))) ' ' units{var_num}],...
    'FontSize',10,'HorizontalAlignment','center');
% clean up
clear SOCAT_grid OAI_grid area_weights t time
% save figure
exportgraphics(gcf,['Figures/' region{n} '_time_series_' var_type{var_num} '.png']);
close