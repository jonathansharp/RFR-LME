
function calculate_oa(vrs,region,num_groups,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

for n = 4%1:length(region)

    % display status
    disp(['Calculating CO2 System (' region{n} ')']);

    % load gridded fCO2, predictors, and models
    load(['Data/LME_Data/' vrs '_prediction_arrays_' region{n}],'LME_prediction');
    load(['Data/RFR-LME/' vrs '_' region{n}],'RFR_LME');
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    load(['Models/' region{n} '/us_lme_models_' vrs '_c' num2str(num_groups(n))],'rfr');

    %% predict TA (and nutrients) using ESPER
    % process ESPER predictors
    lon = repmat(RFR_LME.lon,1,RFR_LME.dim.y,RFR_LME.dim.z);
    lat = repmat(RFR_LME.lat',RFR_LME.dim.x,1,RFR_LME.dim.z);
    lon = lon(LME_prediction.idx);
    lat = lat(LME_prediction.idx);
    depth = ones(size(lon));
    sal = RFR_LME.SSS(LME_prediction.idx);
    tmp = RFR_LME.SST(LME_prediction.idx);
    % predict TA using ESPER-Mixed
    %idx = ~isnan(sal) & ~isnan(tmp);
    [data,u_data] = ESPER_Mixed([1 4 6],[lon lat depth],...
        [sal tmp],[1 2],'Equations',8);
    % pre-allocate TA and nutrients
    RFR_LME.TA = nan(size(LME_prediction.idx));
    RFR_LME.uTA = nan(size(LME_prediction.idx));
    RFR_LME.P = nan(size(LME_prediction.idx));
    RFR_LME.uP = nan(size(LME_prediction.idx));
    RFR_LME.Si = nan(size(LME_prediction.idx));
    RFR_LME.uSi = nan(size(LME_prediction.idx));
    % add TA and nutrients to OA grid
    RFR_LME.TA(LME_prediction.idx) = data.TA;
    RFR_LME.uTA(LME_prediction.idx) = u_data.TA;
    RFR_LME.uP(LME_prediction.idx) = data.phosphate;
    RFR_LME.P(LME_prediction.idx) = u_data.phosphate;
    RFR_LME.Si(LME_prediction.idx) = data.silicate;
    RFR_LME.uSi(LME_prediction.idx) = u_data.silicate;

    %% scale uncertainty over space
    % take gridded absolute delta values
    RFR_LME.ufCO2 = mean(rfr.delta_grid_abs,3,'omitnan');

    % low-pass filter spatial delta fCO2 data twice
%     nan_spc = 1;
    num_cells = 2;
%     while nan_spc > 0 % filter until all grid cells are filled
    for num_filt = 1:2
        RFR_LME.ufCO2 = ... % filter
            smooth2a(RFR_LME.ufCO2,num_cells,num_cells);
%         % check to see if all ocean cells are filled
%         nan_chk = LME_prediction.idx(:,:,1) - ~isnan(RFR_LME.ufCO2);
%         nan_spc = sum(nan_chk(:));
%         num_cells = num_cells+1;
    end
    % filter a third time for PI region
    % RF error was going toward infinity for a few data points
    if n==11
        RFR_LME.ufCO2 = ... % filter
            smooth2a(RFR_LME.ufCO2,num_cells,num_cells);
    end
    % fill unfilled cells with nearest neighbors
    idx = isnan(RFR_LME.ufCO2);
    [lon_grid,lat_grid] = ndgrid(RFR_LME.lon,RFR_LME.lat);
    interp = scatteredInterpolant(lon_grid(~idx),lat_grid(~idx),RFR_LME.ufCO2(~idx),'nearest');
    RFR_LME.ufCO2(idx) = interp(lon_grid(idx),lat_grid(idx));
    RFR_LME.ufCO2(~LME_prediction.idx(:,:,1)) = NaN; % blank out non-ocean cells
    
    % replicate uncertainties over time
    RFR_LME.ufCO2 = ...
        repmat(RFR_LME.ufCO2,1,1,RFR_LME.dim.z);

    % blank ice-filled cells
    RFR_LME.ufCO2(~RFR_LME.idxspc) = NaN;

%     % test plot of spatially scaled absolute delta values
%     figure;
%     h=pcolor(RFR_LME.lon,RFR_LME.lat,...
%         RFR_LME.ufCO2(:,:,1)');
%     set(h,'EdgeColor','none'); colorbar;

    %% scale uncertainty over time
    % load gridded fco2 values
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    % determine annual scaling factors (3-yr to 5-yr periods)
    ann_obs = nan(length(unique(LME.year)),1);
    ann_tot = nan(length(unique(LME.year)),1);
    for y = 1:length(unique(LME.year))
        if y == 1
            ann_obs(y) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,1:36)))));
            ann_tot(y) = sum(sum(sum(LME.idxspc(:,:,1:36))));
        elseif y == 2
            ann_obs(y) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,1:48)))));
            ann_tot(y) = sum(sum(sum(LME.idxspc(:,:,1:48))));
        elseif y == length(unique(LME.year)) - 1
            ann_obs(y) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+48)))));
            ann_tot(y) = sum(sum(sum(LME.idxspc(:,:,(y-3)*12+1:(y-3)*12+48))));
        elseif y == length(unique(LME.year))
            ann_obs(y) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+36)))));
            ann_tot(y) = sum(sum(sum(LME.idxspc(:,:,(y-3)*12+1:(y-3)*12+36))));
        else
            ann_obs(y) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,(y-3)*12+1:(y-3)*12+60)))));
            ann_tot(y) = sum(sum(sum(LME.idxspc(:,:,(y-3)*12+1:(y-3)*12+60))));
        end
            
    end
    ann_per = 100.*(ann_obs./ann_tot);
    % ann_per = repmat(100,length(unique(LME.year)),1);
    % determine seasonal scaling factors (3-month periods)
    mnth_obs = nan(length(unique(LME.month_of_year)),1);
    mnth_tot = nan(length(unique(LME.month_of_year)),1);
    for m = 1:length(unique(LME.month_of_year))
        if m == 1
            mnth_obs(m) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m+11:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(LME.idxspc(:,:,m+11:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m+1:12:end))));
        elseif m == 12
            mnth_obs(m) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m-11:12:end)))));
            mnth_tot(m) = sum(sum(sum(LME.idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m-11:12:end))));
        else
            mnth_obs(m) = sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m-1:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m:12:end)))))+...
                sum(sum(sum(~isnan(LME.fco2_ave_wtd(:,:,m+1:12:end)))));
            mnth_tot(m) = sum(sum(sum(LME.idxspc(:,:,m-1:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m:12:end))))+...
                sum(sum(sum(LME.idxspc(:,:,m+1:12:end))));
        end
    end
    mnth_per = 100.*(mnth_obs./mnth_tot);
    % scale uncertainties
    scaler = nan(LME.dim.z,1);
    for y = 1:length(unique(LME.year))
        for m = 1:length(unique(LME.month_of_year))
            ann_scl = (mean(ann_per)./ann_per(y));
            if ann_scl > 3; ann_scl = 3; end
            mnth_scl = (mean(mnth_per)./mnth_per(m));
            if mnth_scl > 3; mnth_scl = 3; end
            RFR_LME.ufCO2(:,:,(y-1)*12+m) = ...
                RFR_LME.ufCO2(:,:,(y-1)*12+m).*ann_scl.*mnth_scl;
            scaler((y-1)*12+m) = mean([ann_scl mnth_scl]);
        end
    end

    % plot scaling factors
    figure('visible','on');
    plot(datenum(RFR_LME.year,...
        RFR_LME.month_of_year,15),scaler);
    ylabel('Uncertainty Scaler');
    datetick('x');
    exportgraphics(gcf,['Figures/err_scalers_' region{n} '.png']);
    close

    % clean up
    clear SOCAT_grid

    %% calculate other OA Indicators
    fCO2 = RFR_LME.fCO2(LME_prediction.idx);
    carb_system = CO2SYS(data.TA,fCO2,1,5,sal,tmp,NaN,1,NaN,...
        data.silicate,data.phosphate,0,0,1,10,1,2,2);
    % convert -999 to NaN
    carb_system(carb_system==-999) = NaN;
    % pre-allocate OA indicators
    RFR_LME.DIC = nan(size(LME_prediction.idx));
    RFR_LME.pH = nan(size(LME_prediction.idx));
    RFR_LME.OmA = nan(size(LME_prediction.idx));
    RFR_LME.OmC = nan(size(LME_prediction.idx));
    RFR_LME.H = nan(size(LME_prediction.idx));
    RFR_LME.CO3 = nan(size(LME_prediction.idx));
    RFR_LME.RF = nan(size(LME_prediction.idx));
    % add OA indicators to OA grid
    RFR_LME.DIC(LME_prediction.idx) = carb_system(:,2);
    RFR_LME.pH(LME_prediction.idx) = carb_system(:,3);
    RFR_LME.OmA(LME_prediction.idx) = carb_system(:,18);
    RFR_LME.OmC(LME_prediction.idx) = carb_system(:,17);
    RFR_LME.H(LME_prediction.idx) = 10.^-carb_system(:,3);
    RFR_LME.CO3(LME_prediction.idx) = carb_system(:,7);
    RFR_LME.RF(LME_prediction.idx) = carb_system(:,16);
    % clean up
    clear carb_system

    %% calculate other OA Indicator uncertainties
    fCO2_err = RFR_LME.ufCO2(LME_prediction.idx);
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
    RFR_LME.uDIC = nan(size(LME_prediction.idx));
    RFR_LME.upCO2 = nan(size(LME_prediction.idx));
    RFR_LME.upH = nan(size(LME_prediction.idx));
    RFR_LME.uOmA = nan(size(LME_prediction.idx));
    RFR_LME.uOmC = nan(size(LME_prediction.idx));
    RFR_LME.uH = nan(size(LME_prediction.idx));
    RFR_LME.uCO3 = nan(size(LME_prediction.idx));
    RFR_LME.uRF = nan(size(LME_prediction.idx));
    % add OA indicator uncertainties to OA grid
    RFR_LME.uDIC(LME_prediction.idx) = u_carb_system(:,2);
    RFR_LME.upCO2(LME_prediction.idx) = u_carb_system(:,4);
    RFR_LME.upH(LME_prediction.idx) = ... % pH ucertainty by adjusting pH by u[H+]
        RFR_LME.pH(LME_prediction.idx) + ...
        log10(10.^-RFR_LME.pH(LME_prediction.idx) + ...
        (u_carb_system(:,3)./10^9));
    RFR_LME.uOmA(LME_prediction.idx) = u_carb_system(:,11);
    RFR_LME.uOmC(LME_prediction.idx) = u_carb_system(:,10);
    RFR_LME.uH(LME_prediction.idx) = (u_carb_system(:,3)./10^9);
    RFR_LME.uCO3(LME_prediction.idx) = u_carb_system(:,7);
    RFR_LME.uRF(LME_prediction.idx) = u_carb_system(:,9);

    %% save estimated indicators
    if ~isfolder('Data/RFR-LME/'); mkdir('Data/RFR-LME/'); end
    save(['Data/RFR-LME/' vrs '_' region{n}],'RFR_LME','-v7.3');

end

end
