% Import predictor variables

function import_vars(vrs,dpath,source,yr_end,pred_vars_arc)

    % load SOCAT grid
    load(['Data/' vrs '_gridded'],'SOCAT_grid');
    lat = SOCAT_grid.lat;
    lon = SOCAT_grid.lon;
    time = SOCAT_grid.time;
    ocean_mask = SOCAT_grid.percent_sea > 0;
    clear SOCAT_grid

    % 1. obtain sea surface salinity (Options: 'BASS', 'CMEMS')
    import_SSS(dpath,vrs,source.SSS,lat,lon,time,yr_end,'plot_option',0);

    % 2. obtain sea surface height (Options: 'CMEMS', 'NASA')
    import_SSH(dpath,vrs,source.SSH,lat,lon,time,yr_end,'plot_option',0);

    % 3. obtain sea surface temperature
    import_SST(dpath,vrs,source.SST,lat,lon,time,yr_end,'plot_option',0);

    % 4. obtain sea surface ice concentration
    import_IceC(dpath,vrs,source.IceC,lat,lon,time,yr_end,ocean_mask,'plot_option',0);

    % 5. obtain sea surface chlorophyll (Options: 'NASA', 'CMEMS')
    import_CHL(dpath,vrs,source.CHL,lat,lon,time,yr_end,ocean_mask,'plot_option',0);
    
    % 6. Obtain wind speed from ERA5 re-analysis
    import_Wind(dpath,vrs,source.Wind,lat,lon,time,yr_end,'plot_option',0);

    % 7. Obtain mixed layer depth from CMEMS
    import_MLD(dpath,vrs,source.MLD,lat,lon,time,yr_end,'plot_option',0);

    % 8. Obtain atmospheric pressure from NCEP
    import_MSLP(dpath,vrs,source.MSLP,lat,lon,time,yr_end,'plot_option',0);

    % 9. obtain atmospheric pCO2
    import_apCO2(dpath,vrs,source.apCO2,lat,lon,time,yr_end,...
        source.MSLP,source.SSS,source.SST,'plot_option',0);

    % 10. obtain bathymetry from ETOPO2
    import_Bathy(dpath,vrs,source.Bathy,lat,lon,'plot_option',0);

    % Refine ocean mask based on predictor availability
    load(['Data/' vrs '_gridded'],'SOCAT_grid');
    SOCAT_grid.ocean_mask = repmat(ocean_mask,1,1,SOCAT_grid.dim.z);
    for p = 1:length(pred_vars_arc)
        temp_var = ncread(['Data/' pred_vars_arc{p} '_' source.(pred_vars_arc{p}) ...
            '_' vrs '.nc'],pred_vars_arc{p});
        SOCAT_grid.ocean_mask(isnan(temp_var)) = false;
        % sum(sum(SOCAT_grid.ocean_mask(:,:,1),'omitnan'),'omitnan')
        % figure; pcolor(lon,lat,single(SOCAT_grid.ocean_mask(:,:,1))'); shading flat;
    end
    save(['Data/' vrs '_gridded'],'SOCAT_grid','-v7.3');

end