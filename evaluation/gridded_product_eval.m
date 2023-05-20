% Evaluate US LME OA indicators against other gridded products

% Define LMEs
define_regions

%%
for n = 1:length(region)

    %% load CMEMS-LSCE
    CMEMS_LSCE = netcdfreader('/raid/Data/CMEMS_LSCE/CMEMS_LSCE.nc');

    %% remove observations outside general LME limits
    % determine geographic indices
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X)';
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_xmin = find(abs(CMEMS_LSCE.longitude - min(tmp_lon))==...
        min(abs(CMEMS_LSCE.longitude - min(tmp_lon))));
    idx_xmax = find(abs(CMEMS_LSCE.longitude - max(tmp_lon))==...
        min(abs(CMEMS_LSCE.longitude - max(tmp_lon))));
    idx_ymin = find(abs(CMEMS_LSCE.latitude - min(tmp_lat))==...
        min(abs(CMEMS_LSCE.latitude - min(tmp_lat))));
    idx_ymax = find(abs(CMEMS_LSCE.latitude - max(tmp_lat))==...
        min(abs(CMEMS_LSCE.latitude - max(tmp_lat))));
    % pre-allocate
    vars = fields(CMEMS_LSCE);
    % remove gridded observations outside region
    for v = 1:length(vars)
       if ndims(CMEMS_LSCE.(vars{v})) == 3
           CMEMS_LSCE.(region{n}).(vars{v}) = ...
               CMEMS_LSCE.(vars{v})(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
           CMEMS_LSCE = rmfield(CMEMS_LSCE,vars{v});
       end
    end
    % copy other variables
    CMEMS_LSCE.(region{n}).longitude = CMEMS_LSCE.longitude(idx_xmin:idx_xmax);
    CMEMS_LSCE.(region{n}).latitude = CMEMS_LSCE.latitude(idx_ymin:idx_ymax);
    CMEMS_LSCE.(region{n}).time = CMEMS_LSCE.time;
    CMEMS_LSCE = rmfield(CMEMS_LSCE,{'longitude' 'latitude' 'time'});
    % clean up
    clear tmp_lon tmp_lat idx_xmin idx_xmax idx_ymin idx_ymax vars v

    %% remove observations outside refined LME limits
    % determine index based on LME
    CMEMS_LSCE.(region{n}).idxspc = ...
        nan(size(CMEMS_LSCE.(region{n}).spco2,1),size(CMEMS_LSCE.(region{n}).spco2,2));
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    CMEMS_LSCE.(region{n}).idxspc = ...
        inpolygon(...
        repmat(CMEMS_LSCE.(region{n}).longitude,1,length(CMEMS_LSCE.(region{n}).latitude)),...
        repmat(CMEMS_LSCE.(region{n}).latitude',length(CMEMS_LSCE.(region{n}).longitude),1),...
        tmp_lon,tmp_lat);
    CMEMS_LSCE.(region{n}).idxspc = ...
        repmat(CMEMS_LSCE.(region{n}).idxspc,1,1,length(CMEMS_LSCE.(region{n}).time));
    % eliminate gridded data outside LME
    vars = fields(CMEMS_LSCE.(region{n}));
    for v = 1:length(vars)
       if ndims(CMEMS_LSCE.(region{n}).(vars{v})) == 3 && ~strcmp(vars{v},'idxspc')
           CMEMS_LSCE.(region{n}).(vars{v})(~CMEMS_LSCE.(region{n}).idxspc) = NaN;
       end
    end

    %% load US LME data
    load(['Data/' region{n} '/ML_fCO2.mat']);

    %% calculate area-weighted averages
    CMEMS_LSCE



end