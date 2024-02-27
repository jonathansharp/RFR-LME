% Predict TA on grid and perform CO2 system calculations
% 
% This script uses ESPER algorithms to estimate surface TA on a grid for US
% large Marine Ecosystems.
% 
% Written by J.D. Sharp: 1/19/23
% Last updated by J.D. Sharp: 7/27/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions_eiwg

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

    %% scale uncertainty over space
    % take gridded absolute delta values
    OAI_grid.(region{n}).ufCO2 = mean(Val.(region{n}).delta_rfr_grid_abs,3,'omitnan');
    clear Val

    % test plot of original gridded absolute delta values
    % figure;
    % h=pcolor(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,...
    %     OAI_grid.(region{n}).ufCO2');
    % set(h,'EdgeColor','none'); colorbar;
    
    % test with interpolation    
    % OAI_grid.(region{n}).ufCO2_2 = regrid(OAI_grid.(region{n}).lon,...
    %             OAI_grid.(region{n}).lat,OAI_grid.(region{n}).ufCO2,...
    %             OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat);
    % OAI_grid.(region{n}).ufCO2_2(~Preds_grid.(region{n}).idxspc(:,:,1)) = NaN; % blank out non-ocean cells

    % low-pass filter spatial delta fCO2 data twice
%     nan_spc = 1;
    num_cells = 2;
%     while nan_spc > 0 % filter until all grid cells are filled
    for num_filt = 1:2
        OAI_grid.(region{n}).ufCO2 = ... % filter
            smooth2a(OAI_grid.(region{n}).ufCO2,num_cells,num_cells);
%         % check to see if all ocean cells are filled
%         nan_chk = Preds_grid.(region{n}).idxspc(:,:,1) - ~isnan(OAI_grid.(region{n}).ufCO2);
%         nan_spc = sum(nan_chk(:));
%         num_cells = num_cells+1;
    end
    % filter a third time for PI region
    % RF error was going toward infinity for a few data points
    if n==11
        OAI_grid.(region{n}).ufCO2 = ... % filter
            smooth2a(OAI_grid.(region{n}).ufCO2,num_cells,num_cells);
    end
    % fill unfilled cells with nearest neighbors
    idx = isnan(OAI_grid.(region{n}).ufCO2);
    lon_grid = repmat(OAI_grid.(region{n}).lon,1,OAI_grid.(region{n}).dim.y);
    lat_grid = repmat(OAI_grid.(region{n}).lat',OAI_grid.(region{n}).dim.x,1);
    interp = scatteredInterpolant(lon_grid(~idx),lat_grid(~idx),OAI_grid.(region{n}).ufCO2(~idx),'nearest');
    OAI_grid.(region{n}).ufCO2(idx) = interp(lon_grid(idx),lat_grid(idx));
    OAI_grid.(region{n}).ufCO2(~Preds_grid.(region{n}).idxspc(:,:,1)) = NaN; % blank out non-ocean cells
    
    % replicate uncertainties over time
    OAI_grid.(region{n}).ufCO2 = ...
        repmat(OAI_grid.(region{n}).ufCO2,1,1,OAI_grid.(region{n}).dim.z);

    % blank ice-filled cells
    OAI_grid.(region{n}).ufCO2(~OAI_grid.(region{n}).idxspc) = NaN;

%     % test plot of spatially scaled absolute delta values
%     figure;
%     h=pcolor(OAI_grid.(region{n}).lon,OAI_grid.(region{n}).lat,...
%         OAI_grid.(region{n}).ufCO2(:,:,1)');
%     set(h,'EdgeColor','none'); colorbar;

    %% scale uncertainty over time
    % determine annual scaling factors (3-yr to 5-yr periods)
    ann_obs = nan(length(unique(Preds_grid.(region{n}).year)),1);
    ann_tot = nan(length(unique(Preds_grid.(region{n}).year)),1);
    for y = 1:length(unique(Preds_grid.(region{n}).year))
        if y == 1
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,1:36)))));
            ann_tot(y) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,1:36))));
        elseif y == 2
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,1:48)))));
            ann_tot(y) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,1:48))));
        elseif y == length(unique(Preds_grid.(region{n}).year)) - 1
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+48)))));
            ann_tot(y) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+48))));
        elseif y == length(unique(Preds_grid.(region{n}).year))
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+36)))));
            ann_tot(y) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+36))));
        else
            ann_obs(y) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+60)))));
            ann_tot(y) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,(y-3)*12+1:(y-3)*12+60))));
        end
            
    end
    ann_per = 100.*(ann_obs./ann_tot);
    % ann_per = repmat(100,length(unique(Preds_grid.(region{n}).year)),1);
    % determine seasonal scaling factors (3-month periods)
    mnth_obs = nan(length(unique(Preds_grid.(region{n}).month_of_year)),1);
    mnth_tot = nan(length(unique(Preds_grid.(region{n}).month_of_year)),1);
    for m = 1:length(unique(Preds_grid.(region{n}).month_of_year))
        if m == 1
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+11:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m+11:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m+1:12:end))));
        elseif m == 12
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-11:12:end)))));
            mnth_tot(m) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m-11:12:end))));
        else
            mnth_obs(m) = sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(SOCAT_grid.(region{n}).fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m:12:end))))+...
                sum(sum(sum(OAI_grid.(region{n}).idxspc(:,:,m+1:12:end))));
        end
    end
    mnth_per = 100.*(mnth_obs./mnth_tot);
    % scale uncertainties
    scaler = nan(Preds_grid.(region{n}).dim.z,1);
    for y = 1:length(unique(Preds_grid.(region{n}).year))
        for m = 1:length(unique(Preds_grid.(region{n}).month_of_year))
            ann_scl = (mean(ann_per)./ann_per(y));
            if ann_scl > 3; ann_scl = 3; end
            mnth_scl = (mean(mnth_per)./mnth_per(m));
            if mnth_scl > 3; mnth_scl = 3; end
            OAI_grid.(region{n}).ufCO2(:,:,(y-1)*12+m) = ...
                OAI_grid.(region{n}).ufCO2(:,:,(y-1)*12+m).*ann_scl.*mnth_scl;
            scaler((y-1)*12+m) = mean([ann_scl mnth_scl]);
        end
    end

    % plot scaling factors
    figure('visible','off');
    plot(datenum(OAI_grid.(region{n}).year,...
        OAI_grid.(region{n}).month_of_year,15),scaler);
    ylabel('Uncertainty Scaler');
    datetick('x');
    exportgraphics(gcf,['Figures/err_scalers_' region{n} '.png']);
    close

    % clean up
    clear SOCAT_grid

    %% calculate other OA Indicators
    fCO2 = OAI_grid.(region{n}).fCO2(Preds_grid.(region{n}).idxspc);
    carb_system = CO2SYS(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,1,10,1,2,2);
    % convert -999 to NaN
    carb_system(carb_system==-999) = NaN;
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
    % clean up
    clear carb_system

    %% calculate other OA Indicator uncertaintiess
    fCO2_err = OAI_grid.(region{n}).ufCO2(Preds_grid.(region{n}).idxspc);
    u_carb_system = errors(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,u_data.TA,fCO2_err,0,0,...
        u_data.silicate,u_data.phosphate,0,0,'','',0,1,10,1,2,2);
    % tried to pass single precision through but got messed up RF
    % uncertainties. weird rounding?
%     u_carb_system = errors(single(data.TA),single(fCO2),1,5,single(sal),single(tmp),NaN,1,NaN,...
%         single(data.silicate),single(data.phosphate),0,0,single(u_data.TA),single(fCO2_err),0,0,...
%         single(u_data.silicate),single(u_data.phosphate),0,0,'','',0,1,10,1,2,2);
    % convert -999 to NaN
    u_carb_system(u_carb_system==-999) = NaN;
    % pre-allocate OA indicator uncertainties
    OAI_grid.(region{n}).uDIC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).upCO2 = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).upH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uOmA = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uOmC = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uH = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uCO3 = nan(size(Preds_grid.(region{n}).idxspc));
    OAI_grid.(region{n}).uRF = nan(size(Preds_grid.(region{n}).idxspc));
    % add OA indicator uncertainties to OA grid
    OAI_grid.(region{n}).uDIC(Preds_grid.(region{n}).idxspc) = u_carb_system(:,2);
    OAI_grid.(region{n}).upCO2(Preds_grid.(region{n}).idxspc) = u_carb_system(:,4);
    OAI_grid.(region{n}).upH(Preds_grid.(region{n}).idxspc) = ... % pH ucertainty by adjusting pH by u[H+]
        OAI_grid.(region{n}).pH(Preds_grid.(region{n}).idxspc) + ...
        log10(10.^-OAI_grid.(region{n}).pH(Preds_grid.(region{n}).idxspc) + ...
        (u_carb_system(:,3)./10^9));
    OAI_grid.(region{n}).uOmA(Preds_grid.(region{n}).idxspc) = u_carb_system(:,11);
    OAI_grid.(region{n}).uOmC(Preds_grid.(region{n}).idxspc) = u_carb_system(:,10);
    OAI_grid.(region{n}).uH(Preds_grid.(region{n}).idxspc) = (u_carb_system(:,3)./10^9);
    OAI_grid.(region{n}).uCO3(Preds_grid.(region{n}).idxspc) = u_carb_system(:,7);
    OAI_grid.(region{n}).uRF(Preds_grid.(region{n}).idxspc) = u_carb_system(:,9);
    
    %% save estimated OA grid
    save(['Data/' region{n} '/ML_fCO2'],'OAI_grid','-v7.3');

    %% clean up
    clear data carb_system u_carb_system depth fCO2 fCO2_err lat lon
    clear OAI_grid Preds_grid sal TA tmp uTA

end
