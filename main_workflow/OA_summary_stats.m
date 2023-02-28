% Save OA indicator summary statistics
% 
% Written by J.D. Sharp: 2/1/23
% Last updated by J.D. Sharp: 2/23/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' ...
           'HI' 'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};
% variable information
y_span = [1700,2400;250,550;1800,2500;7.9,8.2;1,5;1.5,7.5;6,16;50,300;8,18];
var_type = {'DIC' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 3 2 2 1 1 2];

% preallocate table
stats = nan(length(region),length(var_type)*3);
monthly = nan(288,length(region)*length(var_type));
u_monthly = nan(288,length(region)*length(var_type));
ann = nan(288/12,length(region)*length(var_type));
u_ann = nan(288/12,length(region)*length(var_type));

% loop through each region
for n = 1:length(region)

    for var_num = 1:9

        % load estimated OA grid
        load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
        load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    
        % calculate area-weighted time series
        OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).u_var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
        for t = 1:OAI_grid.(region{n}).dim.z
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

        % calculate annual means
        time_ann = nan(length(OAI_grid.(region{n}).year)/12,1);
        OAI_grid.(region{n}).var_dom_mean_ann = nan(length(OAI_grid.(region{n}).year)/12,1);
        OAI_grid.(region{n}).u_var_dom_mean_ann = nan(length(OAI_grid.(region{n}).year)/12,1);
        for y = 1:length(OAI_grid.(region{n}).year)/12
            time_ann(y) = ...
                mean(time((y-1)*12+1:(y-1)*12+12));
            OAI_grid.(region{n}).var_dom_mean_ann(y) = ...
                mean(OAI_grid.(region{n}).var_dom_mean((y-1)*12+1:(y-1)*12+12));
            OAI_grid.(region{n}).u_var_dom_mean_ann(y) = ...
                mean(OAI_grid.(region{n}).u_var_dom_mean((y-1)*12+1:(y-1)*12+12));
        end
    
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
    
        % fill table with summary statistics
        stats(n,(var_num-1)*3+1) = mean(OAI_grid.(region{n}).var_dom_mean);
        stats(n,(var_num-1)*3+2) = x(2)*12;
        stats(n,(var_num-1)*3+3) = amp;

        % fill table with monthly values
        monthly(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).var_dom_mean;
        u_monthly(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).u_var_dom_mean;

        % fill table with annual values
        ann(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).var_dom_mean_ann;
        u_ann(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).u_var_dom_mean_ann;

        %% plot time series
        % initialize figure
        figure('Visible','off'); hold on;
        set(gcf,'position',[10 10 900 400]);
        % plot
        plot(time,OAI_grid.(region{n}).var_dom_mean,'k-','linewidth',1);
        plot(time_ann,OAI_grid.(region{n}).var_dom_mean_ann,'r-','linewidth',3);
        fill([time;flipud(time)],[OAI_grid.(region{n}).var_dom_mean+OAI_grid.(region{n}).u_var_dom_mean;...
            flipud(OAI_grid.(region{n}).var_dom_mean-OAI_grid.(region{n}).u_var_dom_mean)],'k',...
            'FaceAlpha',0.2,'LineStyle','none');
%         fill([time_ann;flipud(time_ann)],[OAI_grid.(region{n}).var_dom_mean_ann+OAI_grid.(region{n}).u_var_dom_mean_ann;...
%             flipud(OAI_grid.(region{n}).var_dom_mean_ann-OAI_grid.(region{n}).u_var_dom_mean_ann)],'r',...
%             'FaceAlpha',0.2,'LineStyle','none');
        scatter(time_ann,OAI_grid.(region{n}).var_dom_mean_ann,200,'k.');
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
        xL = xlim; yL = ylim;
        % Add recent trend indicator
        if pos_signif_5
            text(xL(2),yL(2),'up','HorizontalAlignment','right','VerticalAlignment','top');
        elseif neg_signif_5
            text(xL(2),yL(2),'down','HorizontalAlignment','right','VerticalAlignment','top');
        else
            text(xL(2),yL(2),'steady','HorizontalAlignment','right','VerticalAlignment','top');
        end
        % Add recent mean indicator
        if pos_mean
            text(xL(2),yL(1),'+','HorizontalAlignment','right','VerticalAlignment','bottom');
        elseif neg_mean
            text(xL(2),yL(1),'-','HorizontalAlignment','right','VerticalAlignment','bottom');
        else
            text(xL(2),yL(1),'o','HorizontalAlignment','right','VerticalAlignment','bottom');
        end
        % save figure
        exportgraphics(gcf,['Figures/' region{n} '_time_series_' var_type{var_num} '.png']);
        close

        % clean up
        clear SOCAT_grid OAI_grid area_weights t time

    end

end

% pre-allocate variable names
VarNameSt = cell(length(var_type)*3,1);
for var_num = 1:9
    VarNameSt{(var_num-1)*3+1} = [var_type{var_num} ',Avg.'];
    VarNameSt{(var_num-1)*3+2} = [var_type{var_num} ',Tr.'];
    VarNameSt{(var_num-1)*3+3} = [var_type{var_num} ',Amp.'];
end
VarNameMnAn = cell(length(region)*length(var_type),1);
for n = 1:length(region)
    for var_num = 1:9
        VarNameMnAn{(n-1)*9+var_num} = [region{n} ',' var_type{var_num}];
    end
end

% convert matrices to tables
stats = array2table(stats,'RowNames',region,'VariableNames',VarNameSt);
monthly = array2table(monthly,'RowNames',arrayfun(@num2str,1:288,'UniformOutput',0),'VariableNames',VarNameMnAn);
ann = array2table(ann,'RowNames',arrayfun(@num2str,1997 + (1:288/12),'UniformOutput',0),'VariableNames',VarNameMnAn);
u_monthly = array2table(u_monthly,'RowNames',arrayfun(@num2str,1:288,'UniformOutput',0),'VariableNames',VarNameMnAn);
u_ann = array2table(u_ann,'RowNames',arrayfun(@num2str,1997 + (1:288/12),'UniformOutput',0),'VariableNames',VarNameMnAn);

% save table
writetable(stats,['Data/SummaryStatistics-' date '.xls'],'WriteRowNames',true);
writetable(monthly,['Data/Monthly-Inds-' date '.xls'],'WriteRowNames',true);
writetable(u_monthly,['Data/Monthly-Inds-Uncer-' date '.xls'],'WriteRowNames',true);
writetable(ann,['Data/Annual-Inds-' date '.xls'],'WriteRowNames',true);
writetable(u_ann,['Data/Annual-Inds-Uncer-' date '.xls'],'WriteRowNames',true);
