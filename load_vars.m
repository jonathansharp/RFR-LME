% Load predictor variables

function load_vars(vrs,dpath)
    % this script defines the bounds of the eleven LMEs
    [lme_shape,lme_idx,region] = define_lme();

    % load SOCAT grid
    load(['Data/' vrs '_gridded'],'SOCAT_grid');

    % display status
    disp(['Loading predictor variables (' region{n} ')']);

    %% Duplicate SOCAT grid for predictors grid
    vars = {'lon' 'lat' 'lim' 'dim' 'month' 'year' 'time' 'month_of_year'};
    for k = 1:length(vars); Preds_grid.(vars{k}) = SOCAT_grid.(vars{k}); end
    clear SOCAT_grid

    %% Calculate distance from coast
    Preds_grid.Dist = ...
        dist2coast(repmat(Preds_grid.lat',Preds_grid.dim.x,1),...
        repmat(Preds_grid.lon,1,Preds_grid.dim.y));

    %% Obtain sea surface salinity from GLORYS
    Preds_grid.SSS = import_SSS(dpath,'BASS',Preds_grid.lat,Preds_grid.lon,Preds_grid.time);

    %% Obtain sea surface height from CMEMS satellite product
%     if n == 5 || n == 6
%         % exclude Bering-Chukchi and Beaufort Seas
%         Preds_grid.(region{n}).SSH = ...
%             nan(size(Preds_grid.(region{n}).idxspc));
%     else
%         Preds_grid.(region{n}).SSH = ...
%             import_SSH_CMEMS(Preds_grid.(region{n}).lat,...
%                             Preds_grid.(region{n}).lon,...
%                             Preds_grid.(region{n}).idxspc(:,:,1),filepath);
%     end

    %% Obtain sea surface temperature from OISSTv2


    %% Obtain sea surface ice concentration from OISSTv2


    %% Obtain sea surface chlorophyll from satellite measurements
%     if n == 5 || n == 6
%         % exclude Bering-Chukchi and Beaufort Seas
%         Preds_grid.(region{n}).CHL = nan(size(Preds_grid.(region{n}).idxspc));
%     else
%         Preds_grid.(region{n}).CHL = ...
%             import_CHL_NASA(Preds_grid.(region{n}).lat,...
%                             Preds_grid.(region{n}).lon,...
%                             Preds_grid.(region{n}).month,...
%                             Preds_grid.(region{n}).idxspc(:,:,1),filepath);
%         Preds_grid.(region{n}).CHL(Preds_grid.(region{n}).CHL<0) = 0.0001;
%     end
    
    %% Obtain wind speed from ERA5 re-analysis

    
    %% Obtain bathymetry from ETOPO2
%     % import_Bathy_ETOPOv2022
%     Preds_grid.(region{n}).Bathy = ...
%         import_Bathy_ETOPOv2022(Preds_grid.(region{n}).lat,...
%                         Preds_grid.(region{n}).lon,...
%                         Preds_grid.(region{n}).idxspc(:,:,1),filepath);
%     % negative trap for bathymetry
%     Preds_grid.(region{n}).Bathy(Preds_grid.(region{n}).Bathy < 0) = 0;

    %% Obtain mixed layer depth

    
    %% Obtain atmospheric pressure from NCEP


    %% Obtain atmospheric pCO2 from NOAA MBL product


    %% Save gridded predictor data
    % save(['Data/' region{n} '/gridded_predictors'],'Preds_grid','-v7.3');

end