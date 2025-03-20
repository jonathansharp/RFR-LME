% Import predictor variables

function import_vars(vrs,dpath,pred_vars,source)

    % load SOCAT grid
    load(['Data/' vrs '_gridded'],'SOCAT_grid');
    lat = SOCAT_grid.lat;
    lon = SOCAT_grid.lon;
    time = SOCAT_grid.time;
    clear SOCAT_grid

    % calculate distance from coast
%     Preds_grid.Dist = ...
%         dist2coast(repmat(Preds_grid.lat',Preds_grid.dim.x,1),...
%         repmat(Preds_grid.lon,1,Preds_grid.dim.y));

    % 1. obtain sea surface salinity (Options: 'BASS', 'GLORYS')
%     import_SSS(dpath,vrs,source{1},lat,lon,time,'plot_option',0);

    % 2. obtain sea surface height (Options: 'CMEMS', 'NASA')
%     import_SSH(dpath,vrs,source{2},lat,lon,time,'plot_option',0);

    % 3. obtain sea surface temperature
%     import_SST(dpath,vrs,source{3},lat,lon,time,'plot_option',0);

    % 4. obtain sea surface ice concentration
    import_IceC(dpath,vrs,source{4},lat,lon,time,'plot_option',0);

    % 5. obtain sea surface chlorophyll (Options: 'NASA', 'CMEMS')
%     import_CHL(dpath,vrs,source{5},lat,lon,time,'plot_option',0);
    
    % Obtain wind speed from ERA5 re-analysis
%     import_SST(dpath,vrs,'CMEMS',lat,lon,time,'plot_option',0);

     % Obtain bathymetry from ETOPO2
%     % import_Bathy_ETOPOv2022
%     Preds_grid.(region{n}).Bathy = ...
%         import_Bathy_ETOPOv2022(Preds_grid.(region{n}).lat,...
%                         Preds_grid.(region{n}).lon,...
%                         Preds_grid.(region{n}).idxspc(:,:,1),filepath);
%     % negative trap for bathymetry
%     Preds_grid.(region{n}).Bathy(Preds_grid.(region{n}).Bathy < 0) = 0;

    % Obtain mixed layer depth
%     import_MLD(dpath,vrs,'CMEMS',lat,lon,time,'plot_option',0);
%     
%     % obtain atmospheric pressure
%     import_MSLP(dpath,vrs,'CMEMS',lat,lon,time,'plot_option',0);
% 
%     % obtain atmospheric pCO2
%     import_SST(dpath,vrs,'CMEMS',lat,lon,time,'plot_option',0);
% 

end