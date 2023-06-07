% Plot OA indicator time series in each LME
% 
% Written by J.D. Sharp: 1/17/23
% Last updated by J.D. Sharp: 2/1/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'PI'};
% variable information
y_span = [1800,2300;250,550;1800,2500;7.9,8.2;1,5;1.5,7.5;6,16;50,300;8,18];
var_type = {'DIC' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 3 2 2 1 1 2];

% loop through variables
for var_num = 1:9

    % initialize figure
    figure('Visible','on'); hold on;
    tcl = tiledlayout(4,3);
    set(gcf,'position',[10 10 1800 1600]);
    nexttile
    set(gca,'Visible','off')

    % add text title
    text(0.5,0.7,['Weighted Mean ' var_lab{var_num} ' Time'],...
        'FontSize',24,'HorizontalAlignment','center');
    text(0.5,0.3,'Series in U.S. LMEs (\muatm)',...
        'FontSize',24,'HorizontalAlignment','center');
%     title(tcl,['Weighted Mean ' var_lab{var_num} ' Time Series in U.S. LMEs (\muatm)'],...
%         'FontSize',24);
    
    % loop through each region
    for n = 1:length(region)
    
        % load estimated OA grid
        load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
        load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    
        % calculate area-weighted time series
        OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).u_var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
        for t = 1:OAI_grid.(region{n}).dim.z
            % remove ice-filled cells
            area_weights(isnan(OAI_grid.(region{n}).(var_type{var_num})(:,:,t))) = NaN;
            OAI_grid.(region{n}).var_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).(var_type{var_num})(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            OAI_grid.(region{n}).u_var_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).(['u' var_type{var_num}])(:,:,t).*...
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
    
        % plot time series
        nexttile
        fill([time;flipud(time)],...
            [OAI_grid.(region{n}).var_dom_mean + OAI_grid.(region{n}).u_var_dom_mean;...
            flipud(OAI_grid.(region{n}).var_dom_mean - OAI_grid.(region{n}).u_var_dom_mean)],...
            rgb('grey'),'LineStyle','none'); hold on
        plot(time,OAI_grid.(region{n}).var_dom_mean,'k','linewidth',2);
        datetick('x');
        xlim([datenum([1998 1 1]) datenum([2022 1 1])]);
        ylim(y_span(var_num,:));
    
        % calculate trend
        [yf,yr,x] = ...
            leastsq2(OAI_grid.(region{n}).month,...
            OAI_grid.(region{n}).var_dom_mean,0,2,[6 12]);
    
        % determine average modelled climatology
        clim = nan(12,1);
        for m = 1:12
            clim(m) = mean(yf(m:12:end));
        end
    
        % calculate amplitude
        amp = max(clim) - min(clim);
    
        % add text to plot
        text(datenum([2010 1 0]),y_span(var_num,2)-(y_span(var_num,2)-y_span(var_num,1))/10,...
            [reg_lab{n} ' | \mu = ' num2str(round(mean(OAI_grid.(region{n}).var_dom_mean),rounder(var_num))) ...
            ' | ' num2str(round(x(2)*12,rounder(var_num))) ...
            ' ' units{var_num} ' yr^{-1} | Amp. = ' ...
            num2str(round(amp,rounder(var_num))) ' ' units{var_num}],...
            'FontSize',10,'HorizontalAlignment','center');
    
        % clean up
        clear SOCAT_grid OAI_grid area_weights t time
    
    end
    
    % save figure
    exportgraphics(gcf,['Figures/all_time_series_' var_type{var_num} '.png']);
    close

end