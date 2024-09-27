% Plot fCO2 ime series in each LME
% 
% Written by J.D. Sharp: 1/17/23
% Last updated by J.D. Sharp: 1/18/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' ...
           'HI' 'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};

% initialize figure
figure('Visible','off'); hold on;
tcl = tiledlayout(3,6);
set(gcf,'position',[10 10 2400 1600]);
title(tcl,'{\itf}CO_{2(RFR)} versus {\itf}CO_{2(meas.)} in U.S. LMEs (\muatm)',...
    'FontSize',24);

% loop through each region
for n = 1:length(region)

    % load estimated fCO2 grid
    load(['Data/' region{n} '/variable_arrays'],'Vars_array');
    load(['Data/' region{n} '/us_lme_model_evals_2'],'Val');

    % plot fCO2 time series
    nexttile
    hold on
    scatter(Vars_array.(region{n}).Y_mod,...
        Val.(region{n}).Y_fit_rfr.all,'.');
    plot([0 1000],[0 1000]);
    xlim([0 1000]);
    ylim([0 1000]);
    hold off

    % add text to plot
    text(50,950,[reg_lab{n} ' | RMSE = ' ...
        num2str(round(Val.(region{n}).rmse_rfr(end),1))],...
        'FontSize',12,'HorizontalAlignment','left');

    % clean up
    clear Val Vars_array

end

% save figure
exportgraphics(gcf,'Figures/all_errs_2.png');
close
