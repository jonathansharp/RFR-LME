% Load predictor variables

function load_vars(vrs,dpath)

    % load SOCAT grid
    load(['Data/' vrs '_gridded'],'SOCAT_grid');

    % duplicate SOCAT grid for predictors grid
    vars = {'lon' 'lat' 'lim' 'dim' 'month' 'year' 'time' 'month_of_year'};
    for k = 1:length(vars); Preds_grid.(vars{k}) = SOCAT_grid.(vars{k}); end
    clear SOCAT_grid

    % calculate distance from coast
    Preds_grid.Dist = ...
        dist2coast(repmat(Preds_grid.lat',Preds_grid.dim.x,1),...
        repmat(Preds_grid.lon,1,Preds_grid.dim.y));

    % obtain sea surface salinity
    Preds_grid.SSS = import_SSS(dpath,vrs,'BASS',Preds_grid.lat,...
        Preds_grid.lon,Preds_grid.time,'plot_option',1);

    % obtain sea surface height
    Preds_grid.SSH = import_SSH(dpath,vrs,'CMEMS',Preds_grid.lat,...
        Preds_grid.lon,Preds_grid.time,'plot_option',1);

    % obtain sea surface temperature
    Preds_grid.SST = import_SST(dpath,vrs,'OISST',Preds_grid.lat,...
        Preds_grid.lon,Preds_grid.time,'plot_option',1);

    % obtain sea surface ice concentration
    Preds_grid.IceC = import_IceC(dpath,vrs,'OISST',Preds_grid.lat,...
        Preds_grid.lon,Preds_grid.time,'plot_option',1);

%     % Obtain sea surface chlorophyll
%     Preds_grid.CHL = import_CHL(dpath,vrs,'CMEMS',Preds_grid.lat,...
%         Preds_grid.lon,Preds_grid.time,'plot_option',0);
%     
%     %% Obtain wind speed from ERA5 re-analysis
%     Preds_grid.Wind = import_SST(dpath,vrs,'CMEMS',Preds_grid.lat,...
%         Preds_grid.lon,Preds_grid.time,'plot_option',0);
% 
%     
%     %% Obtain bathymetry from ETOPO2
% %     % import_Bathy_ETOPOv2022
% %     Preds_grid.(region{n}).Bathy = ...
% %         import_Bathy_ETOPOv2022(Preds_grid.(region{n}).lat,...
% %                         Preds_grid.(region{n}).lon,...
% %                         Preds_grid.(region{n}).idxspc(:,:,1),filepath);
% %     % negative trap for bathymetry
% %     Preds_grid.(region{n}).Bathy(Preds_grid.(region{n}).Bathy < 0) = 0;
% 
%     % obtain mixed layer depth
%     Preds_grid.IceC = import_SST(dpath,vrs,'CMEMS',Preds_grid.lat,...
%         Preds_grid.lon,Preds_grid.time,'plot_option',0);
%     
%     % obtain atmospheric pressure
%     Preds_grid.MSLP = import_MSLP(dpath,vrs,'CMEMS',Preds_grid.lat,...
%         Preds_grid.lon,Preds_grid.time,'plot_option',0);
% 
%     % obtain atmospheric pCO2
%     Preds_grid.pco2 = import_SST(dpath,vrs,'CMEMS',Preds_grid.lat,...
%         Preds_grid.lon,Preds_grid.time,'plot_option',0);
% 
%     % save gridded predictor data
%     % save(['Data/' region{n} '/gridded_predictors'],'Preds_grid','-v7.3');

end