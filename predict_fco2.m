% This function uses random forest regression models to predict fco2 in US
% LMEs and then calculates pCO2 from fCO2.

function predict_fco2(vrs,num_groups,region,varargin)

% process optional inputs
plot_option = 0;
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'plot_option')
        plot_option = varargin{i+1};
    end
end

for n = 1:length(region)

    % display status
    disp(['Predicting fCO2 (' region{n} ')']);
    
    % load variables for prediction    
    load(['Data/LME_Data/' vrs '_' region{n}],'LME');
    load(['Data/LME_Data/' vrs '_GMM_' region{n}],'LME_GMM');
    load(['Data/LME_Data/' vrs '_prediction_arrays_' region{n}],'LME_prediction');

    % load model for prediction
    load(['Models/' region{n} '/us_lme_models_' vrs '_c' num2str(num_groups(n))],'rfr');

    % pre-allocate results
    fco2_rfr_tmp = nan(length(LME_prediction.x),num_groups(n));
    fco2_gmm_probs = nan(length(LME_prediction.x),num_groups(n));

    %% predict fCO2 using RFR for each cluster
    for c = 1:num_groups(n)

        % cluster label
        clab = ['c' num2str(c)];

        % apply random forest regression to data
        temp_rfr = rfr.(clab);
        if strcmp(class(temp_rfr),'TreeBagger')
            fco2_rfr_tmp(:,c) = ...
                predict(temp_rfr,LME_prediction.x);
        else
            fco2_rfr_tmp(:,c) = NaN;
        end

    end

    % average RFR result across clusters
    fco2_rfr = sum(fco2_rfr_tmp.*LME_GMM.probs,2,'omitnan');

    %% re-grid
    % grid information
    RFR_LME.lon = LME.lon;
    RFR_LME.lat = LME.lat;
    RFR_LME.lim = LME.lim;
    RFR_LME.dim = LME.dim;
    RFR_LME.month = LME.month;
    RFR_LME.year = LME.year;
    RFR_LME.month_of_year = LME.month_of_year;
    RFR_LME.idxspc = LME.idxspc;
    RFR_LME.SSS = LME.SSS;
    RFR_LME.SST = LME.SST;

    % assemble fCO2 estimates on grid
    RFR_LME.fCO2 = nan(LME.dim.x,LME.dim.y,LME.dim.z);
    RFR_LME.fCO2(LME_prediction.idx) = fco2_rfr;

    %% process ice-covered grid cells
    % remove fCO2 where ice covers more than 50% of grid cell
    LME_prediction.idx(LME.IceC > 0.5) = 0;
    RFR_LME.fCO2(~LME_prediction.idx) = NaN;

    %% convert fCO2 to pCO2
    % calculate fugacity factor
    TempK = LME.SST + 273.15;
    Delta = 57.7 - 0.118.*TempK;
    b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
    RGasConstant = 83.14462618;
    FugFac = exp((b + 2.*Delta).*LME.MSLP./(RGasConstant.*TempK));
    RFR_LME.pCO2 = RFR_LME.fCO2./FugFac;
    clear TempK Delta b RGasConstant FugFac

    %% plot estimated fCO2 animation
    if plot_option == 1
        create_animation([region{n} '_pCO2'],['RFR-LME_' num2str(num_groups(n))],...
            datenum(RFR_LME.year,RFR_LME.month,15),RFR_LME.lat,RFR_LME.lon,...
            RFR_LME.pCO2,parula,[200 500],'pCO2','(\mumol kg^{-1})');
    end
    
    %% save estimated fCO2 grid
    if ~isfolder('Data/RFR-LME/'); mkdir('Data/RFR-LME/'); end
    save(['Data/RFR-LME/' vrs '_' region{n}],'RFR_LME','-v7.3');

end
