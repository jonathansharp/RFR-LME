% Plot OA indicator climatology in each LME
% 
% Written by J.D. Sharp: 6/7/23
% Last updated by J.D. Sharp: 7/14/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' 'PI'};
% variable information
% y_span = [-100,100;-200,200;-100,100;-0.2,0.2;-1,1;-1,1;6,16;-60,60;-2,2];
var_type = {'DIC' 'fCO2' 'pCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itp}CO_{2}' '{\itf}_{CO2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 1 3 2 2 1 1 2];
co2sys_idx = [2 4 5 1 3 18 17 15 7 16];
deriv_idx = [2 4 5 1 3 11 10 3 7 9];

% loop through variables
for var_num = 1:9

    % preallocate table
    clim_table = nan(length(region),11);

    % initialize figure
    figure('Visible','on'); hold on;
    tcl = tiledlayout(6,2);
    set(gcf,'units','inches','position',[0 0 8.5 11]);

    % add text title
%     text(0.5,0.6,['Weighted Mean ' var_lab{var_num} ' Seasonal'],...
%         'FontSize',12,'HorizontalAlignment','center');
%     text(0.5,0.4,'Cycle in U.S. LMEs (\muatm)',...
%         'FontSize',12,'HorizontalAlignment','center');
%     title(tcl,['Weighted Mean ' var_lab{var_num} ' Time Series in U.S. LMEs (\muatm)'],...
%         'FontSize',24);
    
    % loop through each region
    for n = 1:length(region)
    
        % load estimated OA grid
        load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
        load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
        load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    
        % calculate area-weighted time series
        OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).u_var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).ta_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).dic_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).pco2_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).sst_dom_mean = nan(Preds_grid.(region{n}).dim.z,1);
        OAI_grid.(region{n}).sss_dom_mean = nan(Preds_grid.(region{n}).dim.z,1);
        for t = 1:OAI_grid.(region{n}).dim.z
            area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
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
            OAI_grid.(region{n}).ta_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).TA(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            OAI_grid.(region{n}).dic_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).DIC(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            OAI_grid.(region{n}).pco2_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).pCO2(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            Preds_grid.(region{n}).sst_dom_mean(t) = ...
                squeeze(sum(sum(Preds_grid.(region{n}).SST(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            Preds_grid.(region{n}).sss_dom_mean(t) = ...
                squeeze(sum(sum(Preds_grid.(region{n}).SSS(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            OAI_grid.(region{n}).sil_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).Si(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
            OAI_grid.(region{n}).phos_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).P(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        end
    
        % scale H to nanomoles
        if strcmp(var_type{var_num},'H')
            OAI_grid.(region{n}).var_dom_mean = (10^9).*OAI_grid.(region{n}).var_dom_mean;
            OAI_grid.(region{n}).u_var_dom_mean = (10^9).*OAI_grid.(region{n}).u_var_dom_mean;
        end

        % calculate climatology
        months = (1:12)';
        OAI_grid.(region{n}).var_dom_clim = nan(size(months));
        OAI_grid.(region{n}).u_var_dom_clim = nan(size(months));        
        for m = 1:12
            OAI_grid.(region{n}).var_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).var_dom_mean(m:12:end));
            OAI_grid.(region{n}).u_var_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).u_var_dom_mean(m:12:end));
            OAI_grid.(region{n}).ta_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).ta_dom_mean(m:12:end));
            OAI_grid.(region{n}).dic_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).dic_dom_mean(m:12:end));
            OAI_grid.(region{n}).pco2_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).pco2_dom_mean(m:12:end));
            Preds_grid.(region{n}).sst_dom_clim(m) = ...
                mean(Preds_grid.(region{n}).sst_dom_mean(m:12:end));
            Preds_grid.(region{n}).sss_dom_clim(m) = ...
                mean(Preds_grid.(region{n}).sss_dom_mean(m:12:end));
            OAI_grid.(region{n}).sil_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).sil_dom_mean(m:12:end));
            OAI_grid.(region{n}).phos_dom_clim(m) = ...
                mean(OAI_grid.(region{n}).phos_dom_mean(m:12:end));
        end
    
        % calculate annual means
        var_ann_mean = mean(OAI_grid.(region{n}).var_dom_clim);
        ta_ann_mean = mean(OAI_grid.(region{n}).ta_dom_clim);
        dic_ann_mean = mean(OAI_grid.(region{n}).dic_dom_clim);
        pco2_ann_mean = mean(OAI_grid.(region{n}).pco2_dom_clim);
        sst_ann_mean = mean(Preds_grid.(region{n}).sst_dom_clim);
        sss_ann_mean = mean(Preds_grid.(region{n}).sss_dom_clim);
        sil_ann_mean = mean(OAI_grid.(region{n}).sil_dom_clim);
        phos_ann_mean = mean(OAI_grid.(region{n}).phos_dom_clim);

        % plot time series
        ax=nexttile; hold on
        set(ax,'fontsize',6);
        fill([months;flipud(months)],...
            [(OAI_grid.(region{n}).var_dom_clim-var_ann_mean) + OAI_grid.(region{n}).u_var_dom_clim;...
            flipud((OAI_grid.(region{n}).var_dom_clim-var_ann_mean) - OAI_grid.(region{n}).u_var_dom_clim)],...
            rgb('grey'),'LineStyle','none');
        p1=plot(months,OAI_grid.(region{n}).var_dom_clim-var_ann_mean,'k','linewidth',1);
        xlim([0.5 12.5]); xticks(0:13);
        xticklabels({'' 'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D' ''});
    
        % calculate sensitivity terms
        if strcmp(var_type{var_num},'TA')
            A = derivnum('par1',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,4,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_pco2 = A(:,deriv_idx(var_num));
            A = derivnum('par2',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,4,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_dic = A(:,deriv_idx(var_num));
            A = derivnum('t',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,4,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sst = A(:,deriv_idx(var_num));
            A = derivnum('s',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,4,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sss = A(:,deriv_idx(var_num));
        elseif strcmp(var_type{var_num},'DIC')
            A = derivnum('par1',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).ta_dom_clim,4,1,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_pco2 = A(:,deriv_idx(var_num));
            A = derivnum('par2',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).ta_dom_clim,4,1,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_ta = A(:,deriv_idx(var_num));
            A = derivnum('t',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).ta_dom_clim,4,1,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sst = A(:,deriv_idx(var_num));
            A = derivnum('s',OAI_grid.(region{n}).pco2_dom_clim,...
                OAI_grid.(region{n}).ta_dom_clim,4,1,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sss = A(:,deriv_idx(var_num));
        else
            A = derivnum('par1',OAI_grid.(region{n}).ta_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,1,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_ta = A(:,deriv_idx(var_num));
            A = derivnum('par2',OAI_grid.(region{n}).ta_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,1,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_dic = A(:,deriv_idx(var_num));
            A = derivnum('t',OAI_grid.(region{n}).ta_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,1,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sst = A(:,deriv_idx(var_num));
            A = derivnum('s',OAI_grid.(region{n}).ta_dom_clim,...
                OAI_grid.(region{n}).dic_dom_clim,1,2,Preds_grid.(region{n}).sss_dom_clim,...
                Preds_grid.(region{n}).sst_dom_clim,NaN,0,NaN,OAI_grid.(region{n}).sil_dom_clim,...
                OAI_grid.(region{n}).phos_dom_clim,0,0,1,10,1,2,2);
            d_var_sss = A(:,deriv_idx(var_num));
        end
        
        % adjust derivatives for pH
        if strcmp(var_type{var_num},'pH')
            d_var_ta = OAI_grid.(region{n}).var_dom_clim - ...
                -log10(10.^(-OAI_grid.(region{n}).var_dom_clim)+d_var_ta./(10^9));
            d_var_dic = OAI_grid.(region{n}).var_dom_clim - ...
                -log10(10.^(-OAI_grid.(region{n}).var_dom_clim)+d_var_dic./(10^9));
            d_var_sst = OAI_grid.(region{n}).var_dom_clim - ...
                -log10(10.^(-OAI_grid.(region{n}).var_dom_clim)+d_var_sst./(10^9));
            d_var_sss = OAI_grid.(region{n}).var_dom_clim - ...
                -log10(10.^(-OAI_grid.(region{n}).var_dom_clim)+d_var_sss./(10^9));
        end

        % calculate components from sensitivities
        if strcmp(var_type{var_num},'TA')
            var_dic = d_var_dic.*(OAI_grid.(region{n}).dic_dom_clim'-dic_ann_mean);
            var_sst = d_var_sst.*(Preds_grid.(region{n}).sst_dom_clim'-sst_ann_mean);
            var_sss = d_var_sss.*(Preds_grid.(region{n}).sss_dom_clim'-sss_ann_mean);
            var_pco2 = d_var_pco2.*(OAI_grid.(region{n}).pco2_dom_clim'-pco2_ann_mean);
            var_resid = (OAI_grid.(region{n}).var_dom_clim-var_ann_mean) - ...
                var_pco2 - var_dic - var_sst - var_sss;
        elseif strcmp(var_type{var_num},'DIC')
            var_ta = d_var_ta.*(OAI_grid.(region{n}).ta_dom_clim'-ta_ann_mean);
            var_sst = d_var_sst.*(Preds_grid.(region{n}).sst_dom_clim'-sst_ann_mean);
            var_sss = d_var_sss.*(Preds_grid.(region{n}).sss_dom_clim'-sss_ann_mean);
            var_pco2 = d_var_pco2.*(OAI_grid.(region{n}).pco2_dom_clim'-pco2_ann_mean);
            var_resid = (OAI_grid.(region{n}).var_dom_clim-var_ann_mean) - ...
                var_ta - var_pco2 - var_sst - var_sss;
        else
            var_ta = d_var_ta.*(OAI_grid.(region{n}).ta_dom_clim'-ta_ann_mean);
            var_dic = d_var_dic.*(OAI_grid.(region{n}).dic_dom_clim'-dic_ann_mean);
            var_sst = d_var_sst.*(Preds_grid.(region{n}).sst_dom_clim'-sst_ann_mean);
            var_sss = d_var_sss.*(Preds_grid.(region{n}).sss_dom_clim'-sss_ann_mean);
            var_resid = (OAI_grid.(region{n}).var_dom_clim-var_ann_mean) - ...
                var_ta - var_dic - var_sst - var_sss;
        end


        % plot components
        if strcmp(var_type{var_num},'TA')
            p2=plot(months,var_pco2,'linewidth',1);
            p3=plot(months,var_dic,'linewidth',1);
        elseif strcmp(var_type{var_num},'DIC')
            p2=plot(months,var_ta,'linewidth',1);
            p3=plot(months,var_pco2,'linewidth',1);
        else
            p2=plot(months,var_ta,'linewidth',1);
            p3=plot(months,var_dic,'linewidth',1);
        end
        p4=plot(months,var_sst,'linewidth',1);
        p5=plot(months,var_sss,'linewidth',1);
        p6=plot(months,var_resid,'linewidth',1);
        plot(months,zeros(length(var_ta),1),'k:','linewidth',1);

        % calculate amplitudes
        if strcmp(var_type{var_num},'TA')
            amp_pco2 = max(var_pco2) - min(var_pco2);
            amp_dic = max(var_dic) - min(var_dic);
        elseif strcmp(var_type{var_num},'DIC')
            amp_ta = max(var_ta) - min(var_ta);
            amp_pco2 = max(var_pco2) - min(var_pco2);
        else
            amp_ta = max(var_ta) - min(var_ta);
            amp_dic = max(var_dic) - min(var_dic);
        end
        amp_sst = max(var_sst) - min(var_sst);
        amp_sss = max(var_sss) - min(var_sss);
        amp_resid = max(var_resid) - min(var_resid);
    
        % calculate amplitudes and correlations
        if strcmp(var_type{var_num},'TA')
            r2_dic = corr(var_dic,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sst = corr(var_sst,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sss = corr(var_sss,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_pco2 = corr(var_pco2,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_resid = corr(var_resid,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
        elseif strcmp(var_type{var_num},'DIC')
            r2_ta = corr(var_ta,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sst = corr(var_sst,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sss = corr(var_sss,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_pco2 = corr(var_pco2,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_resid = corr(var_resid,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
        else
            r2_ta = corr(var_ta,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_dic = corr(var_dic,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sst = corr(var_sst,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_sss = corr(var_sss,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
            r2_resid = corr(var_resid,OAI_grid.(region{n}).var_dom_clim-var_ann_mean);
        end

        % add text to plot
        text(6.5,ax.YLim(2)-0.1*(ax.YLim(2)-ax.YLim(1)),...
            reg_lab{n},'FontSize',12,'HorizontalAlignment','center');

        % add legend
        if n == 11
            if strcmp(var_type{var_num},'TA')
                lg=legend([p1 p2 p3 p4 p5 p6],{'Monthly Mean' ...
                    'pCO2 comp.' 'DIC comp.' ...
                    'SST comp.' 'SSS comp.' 'Resid comp.'},...
                    'fontsize',12,'NumColumns',2);
            elseif strcmp(var_type{var_num},'DIC')
                lg=legend([p1 p2 p3 p4 p5 p6],{'Monthly Mean' ...
                    'TA comp.' 'pCO2 comp.' ...
                    'SST comp.' 'SSS comp.' 'Resid comp.'},...
                    'fontsize',12,'NumColumns',2);
            else
                lg=legend([p1 p2 p3 p4 p5 p6],{'Monthly Mean' ...
                    'TA comp.' 'DIC comp.' ...
                    'SST comp.' 'SSS comp.' 'Resid comp.'},...
                    'fontsize',12,'NumColumns',2);
            end
            lg.Position(1) = lg.Position(1)+0.4;
            lg.Position(2) = lg.Position(2)-0.04;
        end

        % add amplitudes and correlations to table
        if strcmp(var_type{var_num},'TA')
            clim_table(n,1) = r2_dic; clim_table(n,2) = amp_dic;
            clim_table(n,3) = r2_pco2; clim_table(n,4) = amp_pco2;
            clim_table(n,5) = r2_sst; clim_table(n,6) = amp_sst;
            clim_table(n,7) = r2_sss; clim_table(n,8) = amp_sss;
            clim_table(n,9) = r2_resid; clim_table(n,10) = amp_resid;
        elseif strcmp(var_type{var_num},'DIC')
            clim_table(n,1) = r2_ta; clim_table(n,2) = amp_ta;
            clim_table(n,3) = r2_pco2; clim_table(n,4) = amp_pco2;
            clim_table(n,5) = r2_sst; clim_table(n,6) = amp_sst;
            clim_table(n,7) = r2_sss; clim_table(n,8) = amp_sss;
            clim_table(n,9) = r2_resid; clim_table(n,10) = amp_resid;
        else
            clim_table(n,1) = r2_ta; clim_table(n,2) = amp_ta;
            clim_table(n,3) = r2_dic; clim_table(n,4) = amp_dic;
            clim_table(n,5) = r2_sst; clim_table(n,6) = amp_sst;
            clim_table(n,7) = r2_sss; clim_table(n,8) = amp_sss;
            clim_table(n,9) = r2_resid; clim_table(n,10) = amp_resid;
        end

        % add true amplitude
        clim_table(n,11) = ...
            max(OAI_grid.(region{n}).var_dom_clim-var_ann_mean) - ...
            min(OAI_grid.(region{n}).var_dom_clim-var_ann_mean);

        % clean up
        clear SOCAT_grid OAI_grid area_weights t months
    
    end
    
    % save figure
    exportgraphics(gcf,['Figures/all_climatology_' var_type{var_num} '.png']);
    close

    % save table
    if strcmp(var_type{var_num},'TA')
        clim_table = array2table(clim_table,'RowNames',region,'VariableNames',...
            {'R2(DIC)' 'Amp.(DIC)' 'R2(pCO2)' 'Amp.(pCO2)' 'R2(SST)' 'Amp.(SST)' ...
            'R2(SSS)' 'Amp.(SSS)' 'R2(resid.)' 'Amp.(resid.)' ' Amp.'});
    elseif strcmp(var_type{var_num},'DIC')
        clim_table = array2table(clim_table,'RowNames',region,'VariableNames',...
            {'R2(TA)' 'Amp.(TA)' 'R2(pCO2)' 'Amp.(pCO2)' 'R2(SST)' 'Amp.(SST)' ...
            'R2(SSS)' 'Amp.(SSS)' 'R2(resid.)' 'Amp.(resid.)' ' Amp.'});
    else
        clim_table = array2table(clim_table,'RowNames',region,'VariableNames',...
            {'R2(TA)' 'Amp.(TA)' 'R2(DIC)' 'Amp.(DIC)' 'R2(SST)' 'Amp.(SST)' ...
            'R2(SSS)' 'Amp.(SSS)' 'R2(resid.)' 'Amp.(resid.)' ' Amp.'});
    end
    writetable(clim_table,['IndsAndStats/ClimatologyStats-' var_type{var_num} '-' date '.xls'],'WriteRowNames',true);

end