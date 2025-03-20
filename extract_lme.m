% Extract LMEs from gridded SOCAT data

function extract_lme(vrs,pred_vars,source,lme_shape,lme_idx,region)

% load SOCAT grid
load(['Data/' vrs '_gridded'],'SOCAT_grid');

% extract each LME from large grid
for n = 1:length(region)

    % display status
    disp(['Extracting ' region{n} ' LME from grid']);

    %% remove observations outside general LME limits
    % determine geographic indices
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X)';
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    idx_xmin = find(abs(SOCAT_grid.lon - min(tmp_lon))==...
        min(abs(SOCAT_grid.lon - min(tmp_lon))));
    idx_xmax = find(abs(SOCAT_grid.lon - max(tmp_lon))==...
        min(abs(SOCAT_grid.lon - max(tmp_lon))));
    idx_ymin = find(abs(SOCAT_grid.lat - min(tmp_lat))==...
        min(abs(SOCAT_grid.lat - min(tmp_lat))));
    idx_ymax = find(abs(SOCAT_grid.lat - max(tmp_lat))==...
        min(abs(SOCAT_grid.lat - max(tmp_lat))));
    % pre-allocate
    vars = fields(SOCAT_grid);
    % remove gridded SOCAT observations outside region
    for v = 1:length(vars)
       if size(SOCAT_grid.(vars{v}),2) == SOCAT_grid.dim.y && ~strcmp(vars{v},'idxspc')
           LME_grid.(region{n}).(vars{v}) = ...
               SOCAT_grid.(vars{v})(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
       end
    end
    % remove gridded predictor observations outside region
    for v = 1:length(pred_vars)
        gridded_predictor = ncread(['Data/' pred_vars{v} '_' ...
            source{v} '_' vrs '.nc'],pred_vars{v});
        LME_grid.(region{n}).(pred_vars{v}) = ...
            gridded_predictor(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
    end
    % copy other variables
    LME_grid.(region{n}).lon = SOCAT_grid.lon(idx_xmin:idx_xmax);
    LME_grid.(region{n}).lat = SOCAT_grid.lat(idx_ymin:idx_ymax);
    LME_grid.(region{n}).month = SOCAT_grid.month;
    LME_grid.(region{n}).year = SOCAT_grid.year;
    LME_grid.(region{n}).month_of_year = SOCAT_grid.month_of_year;
    LME_grid.(region{n}).lim.latmin = min(LME_grid.(region{n}).lat)-0.25;
    LME_grid.(region{n}).lim.latmax = max(LME_grid.(region{n}).lat)+0.25;
    LME_grid.(region{n}).lim.lonmin = min(LME_grid.(region{n}).lon)-0.25;
    LME_grid.(region{n}).lim.lonmax = max(LME_grid.(region{n}).lon)+0.25;
    LME_grid.(region{n}).lim.monthmin = SOCAT_grid.lim.monthmin;
    LME_grid.(region{n}).lim.monthmax = SOCAT_grid.lim.monthmax;
    LME_grid.(region{n}).dim.x = length(LME_grid.(region{n}).lon);
    LME_grid.(region{n}).dim.y = length(LME_grid.(region{n}).lat);
    LME_grid.(region{n}).dim.z = length(LME_grid.(region{n}).month);
    % remove observations outside refined LME limits:
    % determine index based on LME
    LME_grid.(region{n}).idxspc = ...
        nan(LME_grid.(region{n}).dim.x,LME_grid.(region{n}).dim.y);
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X);
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    LME_grid.(region{n}).idxspc = ...
        inpolygon(...
        repmat(LME_grid.(region{n}).lon,1,LME_grid.(region{n}).dim.y),...
        repmat(LME_grid.(region{n}).lat',LME_grid.(region{n}).dim.x,1),...
        tmp_lon,tmp_lat);
    LME_grid.(region{n}).idxspc = ...
        repmat(LME_grid.(region{n}).idxspc,1,1,LME_grid.(region{n}).dim.z);
    % eliminate gridded data outside LME
    vars = fields(LME_grid.(region{n}));
    for v = 1:length(vars)
       if size(LME_grid.(region{n}).(vars{v}),2) == LME_grid.(region{n}).dim.y && ...
               ~strcmp(vars{v},'idxspc')
           if size(LME_grid.(region{n}).(vars{v}),3) == LME_grid.(region{n}).dim.z
               LME_grid.(region{n}).(vars{v})(~LME_grid.(region{n}).idxspc) = NaN;
           else
               LME_grid.(region{n}).(vars{v})(~LME_grid.(region{n}).idxspc(:,:,1)) = NaN;
           end
       end
    end

    % calculate area of LME
    LME_grid.(region{n}).area_tot_km2 = ...
        sum(LME_grid.(region{n}).area_km2(LME_grid.(region{n}).idxspc(:,:,1)).*...
        LME_grid.(region{n}).percent_sea(LME_grid.(region{n}).idxspc(:,:,1)),...
        'omitnan');
    disp(['Area = ' num2str(round(LME_grid.(region{n}).area_tot_km2),0) 'km^2']);

    % clean up
    clear tmp_lon vars v

    % Detrend gridded pCO2 using domain mean:
    % Define area weights
    area_weights = LME_grid.(region{n}).area_km2;
    area_weights = repmat(area_weights,1,1,size(LME_grid.(region{n}).fco2_ave_wtd,3));
    area_weights(isnan(LME_grid.(region{n}).fco2_ave_wtd)) = NaN;
    % Calculate area-weighted domain mean
    LME_grid.(region{n}).fco2_dom_mean = ...
        squeeze(sum(sum(LME_grid.(region{n}).fco2_ave_wtd.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    LME_grid.(region{n}).fco2_dom_mean(LME_grid.(region{n}).fco2_dom_mean == 0) = NaN;
    % Fit trend to area weighted domain mean
    yf = leastsq2(LME_grid.(region{n}).month,...
        LME_grid.(region{n}).fco2_dom_mean,0,0,0);
    % Remove difference from mean for each month
    for m = 1:length(LME_grid.(region{n}).month)
        LME_grid.(region{n}).fco2_ave_wtd_detrend(:,:,m) = ...
            LME_grid.(region{n}).fco2_ave_wtd(:,:,m) + ...
            (mean(yf,'omitnan') - yf(m,:));
    end
    % Calculate area-weighted detrended domain mean
    LME_grid.(region{n}).fco2_dom_mean_detrend = ...
        squeeze(sum(sum(LME_grid.(region{n}).fco2_ave_wtd_detrend.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    LME_grid.(region{n}).fco2_dom_mean_detrend(LME_grid.(region{n}).fco2_dom_mean_detrend == 0) = NaN;

end

% save gridded pco2 and predictor data in individual LMEs
if ~isfolder('Data'); mkdir('Data'); end
save(['Data/' vrs '_gridded_lme'],'LME_grid','-v7.3');

end
