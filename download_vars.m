% Download predictor variables

function download_vars(vrs,dpath,source)

    % load the necessary datasets
    % can this be automated?
    cd([dpath 'CMEMS']);
    path_toolbox = [dpath 'CMEMS/copernicusmarine'];
    datasetID = 'cmems_mod_glo_phy_myint_0.083deg_P1M-m';
    % Create the command to execute    
    command = sprintf('%s subset -i %s -t "%s 12:00:00" -T "%s 12:00:00" -z 0 -Z 0', ...        
        path_toolbox, datasetID, ...
        datestr(datetime(1998,1,1),'yyyy-mm-dd'),...
        datestr(datetime(), 'yyyy-mm-dd'));
    system(command);

end