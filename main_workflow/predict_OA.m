% Predict TA on grid and perform CO2 system calculations
% 
% This script uses ESPER algorithms to estimate surface TA on a grid for US
% large Marine Ecosystems.
% 
% Written by J.D. Sharp: 1/19/23
% Last updated by J.D. Sharp: 4/12/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions

for n = 1:length(region)

    %% display status
    disp(['Calculating CO2 System (' region{n} ')']);

    %% load gridded fCO2, predictors, and models
    load(['Data/' region{n} '/gridded_predictors'],'Preds_grid');
    load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    load(['Data/' region{n} '/us_lme_model_evals'],'Val');

    %% predict TA (and nutrients) using ESPER
    % process ESPER predictors
    lon = repmat(Preds_grid.(region{n}).lon,1,...
        Preds_grid.(region{n}).dim.y,Preds_grid.(region{n}).dim.z);
    lat = repmat(Preds_grid.(region{n}).lat',...
        Preds_grid.(region{n}).dim.x,1,Preds_grid.(region{n}).dim.z);
    lon = lon(Preds_grid.(region{n}).idxspc);
    lat = lat(Preds_grid.(region{n}).idxspc);
    depth = ones(size(lon));
    sal = Preds_grid.(region{n}).SSS(Preds_grid.(region{n}).idxspc);
    tmp = Preds_grid.(region{n}).SST(Preds_grid.(region{n}).idxspc);
    % predict TA using ESPER-Mixed
    [data,u_data] = ESPER_Mixed([1 4 6],[lon lat depth],[sal tmp],[1 2],'Equations',8);
    % pre-allocate TA and nutrients
    OAI_grid.(region{n}).TA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uTA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).P = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uP = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).Si = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uSi = nan(size(Preds_grid.(region{n}).idxspc));
    % add TA and nutrients to OA grid
    OAI_grid.(region{n}).TA(Preds_grid.(region{n}).idxspc) = data.TA;
    OAI_grid.(region{n}).uTA(Preds_grid.(region{n}).idxspc) = u_data.TA;
    OAI_grid.(region{n}).uP(Preds_grid.(region{n}).idxspc) = data.phosphate;
    OAI_grid.(region{n}).P(Preds_grid.(region{n}).idxspc) = u_data.phosphate;
    OAI_grid.(region{n}).Si(Preds_grid.(region{n}).idxspc) = data.silicate;
    OAI_grid.(region{n}).uSi(Preds_grid.(region{n}).idxspc) = u_data.silicate;

    %% plot estimated TA
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).TA,...
        cmocean('haline',24),'TA',...
        'Surface {\itA}_{T} (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).P,...
        cmocean('haline',24),'P',...
        'Surface Phosphate (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).Si,...
        cmocean('haline',24),'Si',...
        'Surface Silicate (\mumol kg^{-1})',region{n});
    % plot TA uncertainties
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
            OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
            OAI_grid.(region{n}).lon,OAI_grid.(region{n}).uTA,...
            cmocean('amp',24),'uTA',...
            'TA Uncertainty (\mumol kg^{-1})',region{n});
    plot_regional_gif(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).uTA,cmocean('amp',20),'uTA',...
        'TA Uncertainty (\mumol kg^{-1})',OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));

    %% scale uncertainty over space (old method)
%     % fit errors to exponential or take mean
%     idx = ~isnan(Val.(region{n}).delta_rfr_grid_abs);
%     Dist_3D = repmat(Preds_grid.(region{n}).Dist,1,1,Preds_grid.(region{n}).dim.z);
%     if n < 11 && n~=3 && n~=4 % exclude AI, EBS, HI, and island regions
%         % fit errors to exponential by distance from coast
%         tbl = table(Dist_3D(idx),Val.(region{n}).delta_rfr_grid_abs(idx));
%         modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x(:,1));
%         beta0 = [1 60 -0.1];
%         mdl = fitnlm(tbl,modelfun,beta0);
%         % plot errors and model
%         figure; hold on;
%         scatter(Dist_3D(idx),Val.(region{n}).delta_rfr_grid_abs(idx));
%         coeffs = mdl.Coefficients{:,'Estimate'};
%         xfit = 0:ceil(max(Dist_3D(idx)));
%         yfit = coeffs(1) + coeffs(2)*exp(coeffs(3).*xfit);
%         plot(xfit,yfit);
%         xlabel('Distance from coast (km)');
%         ylabel('| \Delta{\itf}CO_{2} |');
%         exportgraphics(gcf,['Figures/errs_model_' region{n} '.png']);
%         close
%         % apply exponential to gridded distance
%         OAI_grid.(region{n}).ufCO2 = coeffs(1) + coeffs(2)*exp(coeffs(3).*Dist_3D);
%         OAI_grid.(region{n}).ufCO2(~OAI_grid.(region{n}).idxspc) = NaN; % revove values outside LME
%     else
%         % or use mean error
%         figure; hold on;
%         scatter(Dist_3D(idx),Val.(region{n}).delta_rfr_grid_abs(idx));
%         xfit = 0:ceil(max(Dist_3D(idx)));
%         yfit = repmat(mean(Val.(region{n}).delta_rfr_grid_abs(idx)),length(xfit),1);
%         plot(xfit,yfit);
%         xlabel('Distance from coast (km)');
%         ylabel('| \Delta{\itf}CO_{2} |');
%         exportgraphics(gcf,['Figures/errs_model_' region{n} '.png']);
%         close
%         % or start with mean
%         OAI_grid.(region{n}).ufCO2 = ...
%             repmat(mean(Val.(region{n}).delta_rfr_grid_abs(idx)),size(OAI_grid.(region{n}).fCO2));
%         OAI_grid.(region{n}).ufCO2(~OAI_grid.(region{n}).idxspc) = NaN; % revove values outside LME
%     end

    %% scale uncertainty over space (new method)
    % take gridded absolute delta values
    OAI_grid.(region{n}).ufCO2 = mean(Val.(region{n}).delta_rfr_grid_abs,3,'omitnan');
    % plot of original gridded absolute delta values
    % figure;
    % h=pcolor(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,OAI_grid.(region{n}).ufCO2');
    % set(h,'EdgeColor','none'); colorbar;
    % low-pass filter spatial delta fCO2 data
    nan_spc = 1;
    num_cells = 2;
    while nan_spc > 0 % filter until all grid cells are filled
        OAI_grid.(region{n}).ufCO2 = smooth2a(OAI_grid.(region{n}).ufCO2,num_cells,num_cells);
        OAI_grid.(region{n}).ufCO2(~OAI_grid.(region{n}).idxspc(:,:,1)) = NaN;
        nan_chk = OAI_grid.(region{n}).idxspc(:,:,1) - ~isnan(OAI_grid.(region{n}).ufCO2);
        nan_spc = sum(nan_chk(:));
        num_cells = num_cells+1;
    end
    % plot of spatially scaled absolute delta values
    % figure;
    % h=pcolor(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,OAI_grid.(region{n}).ufCO2');
    % set(h,'EdgeColor','none'); colorbar;
    % replicate uncertainties over time
    OAI_grid.(region{n}).ufCO2 = ...
        repmat(OAI_grid.(region{n}).ufCO2,1,1,OAI_grid.(region{n}).dim.z);

    %% scale uncertainty over time
    % determine annual scaling factors (3-yr to 5-yr periods)
    ann_obs = nan(length(unique(Preds_grid.(region{n}).year)),1);
    ann_tot = nan(length(unique(Preds_grid.(region{n}).year)),1);
    for y = 1:length(unique(Preds_grid.(region{n}).year))
        if y == 1
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,1:36)))));
            ann_tot(y) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,1:36))));
        elseif y == 2
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,1:48)))));
            ann_tot(y) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,1:48))));
        elseif y == length(unique(Preds_grid.(region{n}).year)) - 1
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+48)))));
            ann_tot(y) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+48))));
        elseif y == length(unique(Preds_grid.(region{n}).year))
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+36)))));
            ann_tot(y) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+36))));
        else
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+60)))));
            ann_tot(y) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+60))));
        end
            
    end
    ann_per = 100.*(ann_obs./ann_tot);
    % determine seasonal scaling factors (3-month periods)
    mnth_obs = nan(length(unique(Preds_grid.(region{n}).month_of_year)),1);
    mnth_tot = nan(length(unique(Preds_grid.(region{n}).month_of_year)),1);
    for m = 1:length(unique(Preds_grid.(region{n}).month_of_year))
        if m == 1
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+11:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m+11:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m+1:12:end))));
        elseif m == 12
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-11:12:end)))));
            mnth_tot(m) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m-11:12:end))));
        else
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(SOCAT_grid.(region{n}).idxspc(:,:,m+1:12:end))));
        end
    end
    mnth_per = 100.*(mnth_obs./mnth_tot);
    % scale uncertainties
    scaler = nan(Preds_grid.(region{n}).dim.z,1);
    for y = 1:length(unique(Preds_grid.(region{n}).year))
        for m = 1:length(unique(Preds_grid.(region{n}).month_of_year))
            ann_scl = (mean(ann_per)./ann_per(y));
            if ann_scl > 5; ann_scl = 5; end
            mnth_scl = (mean(mnth_per)./mnth_per(m));
            if mnth_scl > 5; mnth_scl = 5; end
            OAI_grid.(region{n}).ufCO2(:,:,(y-1)*12+m) = ...
                OAI_grid.(region{n}).ufCO2(:,:,(y-1)*12+m).*ann_scl.*mnth_scl;
            scaler((y-1)*12+m) = ann_scl*mnth_scl;
        end
    end
    % plot scaled fCO2 uncertainties
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
            OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
            OAI_grid.(region{n}).lon,OAI_grid.(region{n}).ufCO2,...
            cmocean('amp',24),'ufCO2',...
            '{\itf}CO_{2} Uncertainty (\muatm)',region{n});
    plot_regional_gif(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).ufCO2,cmocean('amp',20),'ufCO2',...
        '{\itf}CO_{2} Uncertainty (\muatm)',OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));

    % plot scaling factors
    figure; plot(1:OAI_grid.(region{n}).dim.z,scaler);
    ylabel('Uncertainty Scaler');
    exportgraphics(gcf,['Figures/err_scalers_' region{n} '.png']);
    close

    %% calculate other OA Indicators
    fCO2 = OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).idxspc);
    fCO2_err = OAI_grid.(region{n}).ufCO2(Preds_grid.(region{n}).idxspc);
    carb_system = CO2SYS(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,1,10,1,2,2);
    u_carb_system = errors(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,u_data.TA,fCO2_err,0,0,...
        u_data.silicate,u_data.phosphate,0,0,'','',0,1,10,1,2,2);
    % convert -999 to NaN
    carb_system(carb_system==-999) = NaN;
    u_carb_system(u_carb_system==-999) = NaN;
    % pre-allocate OA indicators
    OAI_grid.(region{n}).DIC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).pH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).OmA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).OmC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).H = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).CO3 = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).RF = nan(size(Preds_grid.(region{n}).idxspc));
    % add OA indicators to OA grid
    OAI_grid.(region{n}).DIC(Preds_grid.(region{n}).idxspc) = carb_system(:,2);
    OAI_grid.(region{n}).pH(Preds_grid.(region{n}).idxspc) = carb_system(:,3);
    OAI_grid.(region{n}).OmA(Preds_grid.(region{n}).idxspc) = carb_system(:,18);
    OAI_grid.(region{n}).OmC(Preds_grid.(region{n}).idxspc) = carb_system(:,17);
    OAI_grid.(region{n}).H(Preds_grid.(region{n}).idxspc) = 10.^-carb_system(:,3);
    OAI_grid.(region{n}).CO3(Preds_grid.(region{n}).idxspc) = carb_system(:,7);
    OAI_grid.(region{n}).RF(Preds_grid.(region{n}).idxspc) = carb_system(:,16);
    % pre-allocate OA indicator uncertainties
    OAI_grid.(region{n}).uDIC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).upH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uOmA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uOmC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uCO3 = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uRF = nan(size(Preds_grid.(region{n}).idxspc));
    % add OA indicator uncertainties to OA grid
    OAI_grid.(region{n}).uDIC(Preds_grid.(region{n}).idxspc) = u_carb_system(:,2);
    OAI_grid.(region{n}).upH(Preds_grid.(region{n}).idxspc) = ... % pH ucertainty by adjusting pH by u[H+]
        carb_system(:,3) + log10(10.^-carb_system(:,3)+(u_carb_system(:,3)./10^9));
    OAI_grid.(region{n}).uOmA(Preds_grid.(region{n}).idxspc) = u_carb_system(:,11);
    OAI_grid.(region{n}).uOmC(Preds_grid.(region{n}).idxspc) = u_carb_system(:,10);
    OAI_grid.(region{n}).uH(Preds_grid.(region{n}).idxspc) = 10.^-u_carb_system(:,3);
    OAI_grid.(region{n}).uCO3(Preds_grid.(region{n}).idxspc) = u_carb_system(:,7);
    OAI_grid.(region{n}).uRF(Preds_grid.(region{n}).idxspc) = u_carb_system(:,9);
    
    %% plot estimated OA indicators
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).DIC,...
        flipud(jet(20)),'DIC','Dissolved Inorganic Carbon (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).pH,...
        flipud(jet(20)),'pH','Surface pH_{T}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).OmA,...
        flipud(jet(16)),'OmA','Surface \Omega_{A}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).OmC,...
        flipud(jet(16)),'OmC','Surface \Omega_{C}',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,(10^9).*OAI_grid.(region{n}).H,...
        jet(20),'H','Surface [H^{+}] (nmol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).CO3,...
        jet(15),'CO3','Surface [CO_{3}^{2-}] (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).RF,...
        flipud(jet(28)),'RF','Surface RF',region{n});
    plot_regional_gif(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).TA,cmocean('haline',24),'TA',...
        'Surface {\itA}_{T} (\mumol kg^{-1})',OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));
    % plot pH uncertainties
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
            OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
            OAI_grid.(region{n}).lon,OAI_grid.(region{n}).upH,...
            cmocean('amp',20),'upH',...
            'pH Uncertainty',region{n});
    plot_regional_gif(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).lat,OAI_grid.(region{n}).lon,...
        OAI_grid.(region{n}).upH,cmocean('amp',20),'upH',...
        'pH Uncertainty',OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,region{n},lme_shape(lme_idx.(region{n})));

    %% plot estimated OA indicator uncertainties
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).uTA,...
        flipud(jet(28)),'uTA','Surface TA Uncertainty (\mumol kg^{-1})',region{n});
    plot_temporal_mean(OAI_grid.(region{n}).lim,...
        OAI_grid.(region{n}).dim,OAI_grid.(region{n}).lat,...
        OAI_grid.(region{n}).lon,OAI_grid.(region{n}).upH,...
        flipud(jet(28)),'upH','Surface pH Uncertainty',region{n});

    %% save estimated OA grid
    save(['Data/' region{n} '/ML_fCO2'],'OAI_grid','-v7.3');

    %% clean up
    clear carb_system u_carb_system depth fCO2 lat lon OAI_grid Preds_grid sal TA tmp uTA

end

%% plot time series
OA_summary_stats
OA_time_series

%% plot OA indicators across full region
% plot_temporal_mean_full(1900,2500,25,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
% plot_temporal_mean_full(1700,2300,25,flipud(jet(15)),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
% plot_temporal_mean_full(7.9,8.2,0.02,jet(15),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
% plot_temporal_mean_full(0,5,0.25,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
% plot_temporal_mean_full(1,6,0.25,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
% plot_temporal_mean_full(5,10,0.25,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
% plot_temporal_mean_full(50,250,10,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
% plot_temporal_mean_full(8,15,0.5,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)
plot_temporal_mean_full(0,40,2,cmocean('amp',21),'ufCO2','Sea Surface {\itf}CO_{2} Uncertainty',region,lme_shape,lme_idx)

%% plot gifs of OA indicators across full region
% plot_full_gif(1900,2500,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
% plot_full_gif(1700,2300,jet(15),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
% plot_full_gif(7.9,8.2,flipud(jet(15)),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
% plot_full_gif(0,5,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
% plot_full_gif(1,6,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
% plot_full_gif(5,10,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
% plot_full_gif(50,250,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
% plot_full_gif(8,16,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)
plot_full_gif(0,40,cmocean('amp',21),'ufCO2','Sea Surface {\itf}CO_{2} Uncertainty',region,lme_shape,lme_idx)

%% plot OA indicators across full region (seasonally)
% plot_temporal_mean_full_seas(1900,2500,25,cmocean('haline',24),'TA','Sea Surface {\itA}_{T} (\mumol kg^{-1})',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(1700,2300,25,flipud(jet(15)),'DIC','Dissolved Inorganic Carbon',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(7.9,8.2,0.02,flipud(jet(15)),'pH','Sea Surface pH_{T}',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(0,5,0.25,flipud(jet(20)),'OmA','Sea Surface \Omega_{A}',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(1,6,0.25,flipud(jet(20)),'OmC','Sea Surface \Omega_{C}',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(5,10,0.25,parula(20),'H','Sea Surface [H^{+}]',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(50,250,10,parula(20),'CO3','Sea Surface [CO_{3}^{2-}]',region,lme_shape,lme_idx)
% plot_temporal_mean_full_seas(8,15,0.5,parula(14),'RF','Sea Surface RF',region,lme_shape,lme_idx)
