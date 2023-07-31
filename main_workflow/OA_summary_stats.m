% Save OA indicator summary statistics
% 
% Written by J.D. Sharp: 2/1/23
% Last updated by J.D. Sharp: 7/14/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'PI'};
% variable information
var_type = {'DIC' 'pCO2' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF' 'SST' 'SSS' 'TA_DIC'};
units = {'\mumol kg^{-1}' '\muatm' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' '' 'degC' '' ''};
rounder = [1 1 1 1 3 2 2 1 1 2 1 1 1];

% preallocate tables
stats = nan(length(region),length(var_type)*5);
monthly = nan(288,length(region)*length(var_type));
u_monthly = nan(288,length(region)*length(var_type));
ann = nan(288/12,length(region)*length(var_type));
u_ann = nan(288/12,length(region)*length(var_type));
means_table = nan(length(region),length(var_type));
iav_table = nan(length(region),length(var_type));
trends_table = nan(length(region),length(var_type)*2);
amp_table = nan(length(region),length(var_type));

% loop through each region
for n = 1:length(region)

    for var_num = 3%1:length(var_type)

        % load estimated OA grid
        load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
        load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
        load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    
        % calculate area-weighted time series
        OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).u_var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
        for t = 1:OAI_grid.(region{n}).dim.z
            % establish non-ice-covered fraction
            open_per = sum(OAI_grid.(region{n}).idxspc(:,:,t),'all')./...
                sum(SOCAT_grid.(region{n}).idxspc(:,:,t),'all');
            % remove ice-filled cells from area weights
            area_weights(~OAI_grid.(region{n}).idxspc(:,:,t)) = NaN;
            % calculate area-weighted means
            if strcmp(var_type{var_num},'SST') || strcmp(var_type{var_num},'SSS')
                OAI_grid.(region{n}).var_dom_mean(t) = ...
                    squeeze(sum(sum(Preds_grid.(region{n}).(var_type{var_num})(:,:,t).*...
                        area_weights,1,'omitnan'),2,'omitnan'))./...
                        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
                OAI_grid.(region{n}).u_var_dom_mean(t) = NaN;
            elseif strcmp(var_type{var_num},'TA_DIC')
                OAI_grid.(region{n}).var_dom_mean(t) = ...
                    squeeze(sum(sum((OAI_grid.(region{n}).TA(:,:,t)./...
                    OAI_grid.(region{n}).DIC(:,:,t)).*...
                        area_weights,1,'omitnan'),2,'omitnan'))./...
                        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
                OAI_grid.(region{n}).u_var_dom_mean(t) = NaN;
            else
                OAI_grid.(region{n}).var_dom_mean(t) = ...
                    squeeze(sum(sum(OAI_grid.(region{n}).(var_type{var_num})(:,:,t).*...
                        area_weights,1,'omitnan'),2,'omitnan'))./...
                        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
                OAI_grid.(region{n}).u_var_dom_mean(t) = ...
                    squeeze(sum(sum(OAI_grid.(region{n}).(['u' var_type{var_num}])(:,:,t).*...
                        area_weights,1,'omitnan'),2,'omitnan'))./...
                        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));                
            end
            % remove means when region is >50% ice (retaining these for now)
%             if open_per < 0.5
%                 OAI_grid.(region{n}).var_dom_mean(t) = NaN;
%                 OAI_grid.(region{n}).u_var_dom_mean(t) = NaN;
%             end
        end
    
        % scale H and uH fom moles to nanomoles
        if strcmp(var_type{var_num},'H')
            OAI_grid.(region{n}).var_dom_mean = (10^9).*OAI_grid.(region{n}).var_dom_mean;
            OAI_grid.(region{n}).u_var_dom_mean = (10^9).*OAI_grid.(region{n}).u_var_dom_mean;
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
                mean(OAI_grid.(region{n}).var_dom_mean((y-1)*12+1:(y-1)*12+12),'omitnan');
            OAI_grid.(region{n}).u_var_dom_mean_ann(y) = ...
                mean(OAI_grid.(region{n}).u_var_dom_mean((y-1)*12+1:(y-1)*12+12),'omitnan');
        end
    
        % calculate long-term mean
        mean_lt = mean(OAI_grid.(region{n}).var_dom_mean,'omitnan');

        % calculate trend
        [yf,yr,x,err] = ...
            leastsq2(1:length(OAI_grid.(region{n}).year)/12,...
            OAI_grid.(region{n}).var_dom_mean_ann,0,0,0);

        % uncertainty on trend
        [acov,acor,lag,dof] = ...
            autocov2(1:length(OAI_grid.(region{n}).year)/12,yr,10);
        edof = dof - 2; % subtract number of parameters to get effective dof
        tr_uncer = ... % scale uncertainty using edof
            err(2)*(sqrt(length(OAI_grid.(region{n}).year)/12)./sqrt(edof));
    
        % calculate interannual variability
        iav = std(yr);

        % determine average modelled climatology
        clim = nan(12,1);
        clim_uncer = nan(12,1);
        for m = 1:12
            clim(m) = mean(OAI_grid.(region{n}).var_dom_mean(m:12:end),'omitnan');
            clim_uncer(m) = mean(OAI_grid.(region{n}).u_var_dom_mean(m:12:end),'omitnan');
        end
    
        % calculate amplitude
        if sum(~isnan(clim)) > 8
            amp = max(clim) - min(clim);
        else
            amp = NaN;
        end

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
        stats(n,(var_num-1)*5+1) = mean_lt; % mean
        stats(n,(var_num-1)*5+2) = x(2); % trend
        stats(n,(var_num-1)*5+3) = tr_uncer; % trend uncertainty
        stats(n,(var_num-1)*5+4) = amp; % amplitude
        stats(n,(var_num-1)*5+5) = iav; % interannual variability

        % fill individual tables
        means_table(n,var_num) = mean_lt;
        trends_table(n,(var_num-1)*2+1) = x(2);
        trends_table(n,(var_num-1)*2+2) = tr_uncer;
        amp_table(n,var_num) = amp;
        iav_table(n,var_num) = iav;

        % fill table with monthly values
        monthly(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).var_dom_mean;
        u_monthly(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).u_var_dom_mean;

        % fill table with annual values
        ann(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).var_dom_mean_ann;
        u_ann(:,(n-1)*length(var_type)+var_num) = OAI_grid.(region{n}).u_var_dom_mean_ann;

        %% plot time series
        % initialize figure
        figure('Visible','on'); hold on;
        set(gcf,'position',[10 10 900 400]);
        % plot
        plot(time,OAI_grid.(region{n}).var_dom_mean,'k-','linewidth',1);
        plot(time_ann,OAI_grid.(region{n}).var_dom_mean_ann,'r-','linewidth',3);
%         fill([time;flipud(time)],[OAI_grid.(region{n}).var_dom_mean+OAI_grid.(region{n}).u_var_dom_mean;...
%             flipud(OAI_grid.(region{n}).var_dom_mean-OAI_grid.(region{n}).u_var_dom_mean)],'k',...
%             'FaceAlpha',0.2,'LineStyle','none');
        fill([time_ann;flipud(time_ann)],[OAI_grid.(region{n}).var_dom_mean_ann+OAI_grid.(region{n}).u_var_dom_mean_ann;...
            flipud(OAI_grid.(region{n}).var_dom_mean_ann-OAI_grid.(region{n}).u_var_dom_mean_ann)],'k',...
            'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
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
        title([reg_lab{n} ' | \mu = ' num2str(round(mean(OAI_grid.(region{n}).var_dom_mean,'omitnan'),rounder(var_num))) ...
            ' ' units{var_num} ' | Tr. = ' num2str(round(x(2),rounder(var_num))) ...
            ' ' units{var_num} ' yr^{-1} | Amp. = ' ...
            num2str(round(amp,rounder(var_num))) ' ' units{var_num}],...
            'FontSize',10,'HorizontalAlignment','center');
        xL = xlim; yL = ylim;
        % Add recent trend indicator
%         if pos_signif_5
%             text(xL(2),yL(2),'up','HorizontalAlignment','right','VerticalAlignment','top');
%         elseif neg_signif_5
%             text(xL(2),yL(2),'down','HorizontalAlignment','right','VerticalAlignment','top');
%         else
%             text(xL(2),yL(2),'steady','HorizontalAlignment','right','VerticalAlignment','top');
%         end
        % Add recent mean indicator
%         if pos_mean
%             text(xL(2),yL(1),'+','HorizontalAlignment','right','VerticalAlignment','bottom');
%         elseif neg_mean
%             text(xL(2),yL(1),'-','HorizontalAlignment','right','VerticalAlignment','bottom');
%         else
%             text(xL(2),yL(1),'o','HorizontalAlignment','right','VerticalAlignment','bottom');
%         end
        % save figure
        exportgraphics(gcf,['Figures/' region{n} '_time_series_' var_type{var_num} '.png']);
        close

        %% plot climatology
        % initialize figure
        figure('Visible','on'); hold on;
        set(gcf,'position',[10 10 700 400]);
        % plot
        plot(1:12,clim,'k-','linewidth',3);
        fill([1:12,fliplr(1:12)]',[clim+clim_uncer;...
            flipud(clim-clim_uncer)],'k',...
            'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
        scatter(1:12,clim,200,'k.');
        xlim([0.5 12.5]);
        xticks(1:12);
        xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
        % add text to plot
        title([reg_lab{n} ' | \mu = ' num2str(round(mean(OAI_grid.(region{n}).var_dom_mean),rounder(var_num))) ...
            ' ' units{var_num} ' | Amp. = ' ...
            num2str(round(amp,rounder(var_num))) ' ' units{var_num}],...
            'FontSize',10,'HorizontalAlignment','center');
        % Add recent trend indicator
%         if pos_signif_5
%             text(xL(2),yL(2),'up','HorizontalAlignment','right','VerticalAlignment','top');
%         elseif neg_signif_5
%             text(xL(2),yL(2),'down','HorizontalAlignment','right','VerticalAlignment','top');
%         else
%             text(xL(2),yL(2),'steady','HorizontalAlignment','right','VerticalAlignment','top');
%         end
        % Add recent mean indicator
%         if pos_mean
%             text(xL(2),yL(1),'+','HorizontalAlignment','right','VerticalAlignment','bottom');
%         elseif neg_mean
%             text(xL(2),yL(1),'-','HorizontalAlignment','right','VerticalAlignment','bottom');
%         else
%             text(xL(2),yL(1),'o','HorizontalAlignment','right','VerticalAlignment','bottom');
%         end
        % save figure
        exportgraphics(gcf,['Figures/' region{n} '_climatology_' var_type{var_num} '.png']);
        close

        %% save LME-specific indicator time series
        if ~exist(['IndsAndStats/' region{n}],'dir')
            mkdir(['IndsAndStats/' region{n}]);
        end
        mnth_tmp = [OAI_grid.(region{n}).year,OAI_grid.(region{n}).month_of_year,...
            OAI_grid.(region{n}).var_dom_mean,OAI_grid.(region{n}).u_var_dom_mean];
        writetable(array2table(mnth_tmp,'VariableNames',{'Year' 'Month' 'Value' 'Uncertainty'}),...
            ['IndsAndStats/' region{n} '/monthly_' var_type{var_num} '.csv']);
        ann_tmp = [unique(OAI_grid.(region{n}).year),...
            OAI_grid.(region{n}).var_dom_mean_ann,OAI_grid.(region{n}).u_var_dom_mean_ann];
        writetable(array2table(ann_tmp,'VariableNames',{'Year' 'Value' 'Uncertainty'}),...
            ['IndsAndStats/' region{n} '/annual_' var_type{var_num} '.csv'],'WriteRowNames',true);

        %% clean up
        clear SOCAT_grid OAI_grid area_weights t time

        %% save LME-specific indicator time series

    end

end

% pre-allocate variable names
VarNameSt = cell(length(var_type)*5,1);
for var_num = 1:length(var_type)
    VarNameSt{(var_num-1)*5+1} = [var_type{var_num} ',Avg.'];
    VarNameSt{(var_num-1)*5+2} = [var_type{var_num} ',Tr.'];
    VarNameSt{(var_num-1)*5+3} = [var_type{var_num} ',Tr. Uncer.'];
    VarNameSt{(var_num-1)*5+4} = [var_type{var_num} ',Amp.'];
    VarNameSt{(var_num-1)*5+5} = [var_type{var_num} ',IAV'];
end
VarNameMnAn = cell(length(region)*length(var_type),1);
for n = 1:length(region)
    for var_num = 1:length(var_type)
        VarNameMnAn{(n-1)*length(var_type)+var_num} = [region{n} ',' var_type{var_num}];
    end
end

% convert matrices to tables
stats = array2table(stats,'RowNames',region,'VariableNames',VarNameSt);
means_table = array2table(means_table,'RowNames',region,'VariableNames',VarNameSt(1:5:end));
trend_var_names = [VarNameSt(2:5:end),VarNameSt(3:5:end)]';
trends_table = array2table(trends_table,'RowNames',region,'VariableNames',trend_var_names);
amp_table = array2table(amp_table,'RowNames',region,'VariableNames',VarNameSt(4:5:end));
iav_table = array2table(iav_table,'RowNames',region,'VariableNames',VarNameSt(5:5:end));
monthly = array2table(monthly,'RowNames',arrayfun(@num2str,1:288,'UniformOutput',0),'VariableNames',VarNameMnAn);
ann = array2table(ann,'RowNames',arrayfun(@num2str,1997 + (1:288/12),'UniformOutput',0),'VariableNames',VarNameMnAn);
u_monthly = array2table(u_monthly,'RowNames',arrayfun(@num2str,1:288,'UniformOutput',0),'VariableNames',VarNameMnAn);
u_ann = array2table(u_ann,'RowNames',arrayfun(@num2str,1997 + (1:288/12),'UniformOutput',0),'VariableNames',VarNameMnAn);

% save table
writetable(stats,['IndsAndStats/SummaryStatistics-' date '.xls'],'WriteRowNames',true);
writetable(means_table,['IndsAndStats/MeansTable-' date '.xls'],'WriteRowNames',true);
writetable(trends_table,['IndsAndStats/TrendsTable-' date '.xls'],'WriteRowNames',true);
writetable(amp_table,['IndsAndStats/AmplitudeTable-' date '.xls'],'WriteRowNames',true);
writetable(iav_table,['IndsAndStats/IAVTable-' date '.xls'],'WriteRowNames',true);
writetable(monthly,['IndsAndStats/Monthly-Inds-' date '.xls'],'WriteRowNames',true);
writetable(u_monthly,['IndsAndStats/Monthly-Inds-Uncer-' date '.xls'],'WriteRowNames',true);
writetable(ann,['IndsAndStats/Annual-Inds-' date '.xls'],'WriteRowNames',true);
writetable(u_ann,['IndsAndStats/Annual-Inds-Uncer-' date '.xls'],'WriteRowNames',true);
