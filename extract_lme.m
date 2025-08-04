% Extract LMEs from gridded SOCAT data

function extract_lme(vrs,pred_vars,pred_vars_arc,source,lme_shape,lme_idx,region)

% load SOCAT grid
load(['Data/' vrs '_gridded'],'SOCAT_grid');
% SOCAT_vars = fieldnames(SOCAT_grid);
% keep_vars = {'fco2_ave_wtd' 'lat' 'lon' 'dim' 'lim' 'month' ...
%     'year' 'month_of_year' 'area_km2' 'percent_sea'};
% idx = false(size(SOCAT_vars));
% for p = 1:length(keep_vars)
%     idx = idx | strcmp(SOCAT_vars,keep_vars{p});
% end
% SOCAT_grid = ...
%     rmfield(SOCAT_grid,SOCAT_vars(~idx));

% load predictor variables
for p = 1:length(pred_vars)
    SOCAT_grid.(pred_vars{p}) = ...
        ncread(['Data/' pred_vars{p} '_' source.(pred_vars{p}) '_' vrs '.nc'],pred_vars{p});
end

% apply mask to each predictor variable according to bathymetry
% per_sea = repmat(SOCAT_grid.percent_sea,1,1,SOCAT_grid.dim.z);
% ocean_mask = true(SOCAT_grid.dim.x,SOCAT_grid.dim.y,SOCAT_grid.dim.z);
% ocean_mask(per_sea==0) = false;
% for p = 1:length(pred_vars_arc)
%     if ndims(SOCAT_grid.(pred_vars_arc{p})) == 3
%         ocean_mask(isnan(SOCAT_grid.(pred_vars_arc{p}))) = false;
%     end
% end
% for p = 1:length(pred_vars)
%     if ndims(SOCAT_grid.(pred_vars{p})) == 2
%         SOCAT_grid.(pred_vars{p})(~ocean_mask(:,:,1)) = NaN;
%     else
%         SOCAT_grid.(pred_vars{p})(~ocean_mask) = NaN;
%     end
%     sum(sum(isnan(SOCAT_grid.(pred_vars{p})(:,:,1))))
%     % figure; pcolor(SOCAT_grid.lon,SOCAT_grid.lat,SOCAT_grid.SSH(:,:,1)'); shading flat; colorbar;
% end

% extract each LME from large grid
for n = 1:length(region)

    % display status
    disp(['Extracting ' region{n} ' LME from grid']);

    %% remove observations outside general LME limits
    % determine geographic indices
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X,'format','0-360')';
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
           LME.(vars{v}) = ...
               SOCAT_grid.(vars{v})(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
       end
    end
    % remove gridded predictor observations outside region
    for v = 1:length(pred_vars)
        gridded_predictor = ncread(['Data/' pred_vars{v} '_' ...
            source.(pred_vars{v}) '_' vrs '.nc'],pred_vars{v});
        LME.(pred_vars{v}) = ...
            gridded_predictor(idx_xmin:idx_xmax,idx_ymin:idx_ymax,:);
    end
    % copy other variables
    LME.lon = SOCAT_grid.lon(idx_xmin:idx_xmax);
    LME.lat = SOCAT_grid.lat(idx_ymin:idx_ymax);
    LME.month = SOCAT_grid.month;
    LME.year = SOCAT_grid.year;
    LME.month_of_year = SOCAT_grid.month_of_year;
    LME.lim.latmin = min(LME.lat)-0.25;
    LME.lim.latmax = max(LME.lat)+0.25;
    LME.lim.lonmin = min(LME.lon)-0.25;
    LME.lim.lonmax = max(LME.lon)+0.25;
    LME.lim.monthmin = SOCAT_grid.lim.monthmin;
    LME.lim.monthmax = SOCAT_grid.lim.monthmax;
    LME.dim.x = length(LME.lon);
    LME.dim.y = length(LME.lat);
    LME.dim.z = length(LME.month);
    % remove observations outside refined LME limits:
    % determine index based on LME
    LME.idxspc = nan(LME.dim.x,LME.dim.y);
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X,'format','0-360');
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    LME.idxspc = inpolygon(repmat(LME.lon,1,LME.dim.y),...
        repmat(LME.lat',LME.dim.x,1),tmp_lon,tmp_lat);
    LME.idxspc = repmat(LME.idxspc,1,1,LME.dim.z);
    LME.idxspc(~LME.ocean_mask) = false;
    % eliminate gridded data outside LME
    vars = fields(LME);
    for v = 1:length(vars)
       if size(LME.(vars{v}),2) == LME.dim.y && ...
               ~strcmp(vars{v},'idxspc') && ...
               ~strcmp(vars{v},'ocean_mask')
           if size(LME.(vars{v}),3) == LME.dim.z
               LME.(vars{v})(~LME.idxspc) = NaN;
           elseif size(LME.(vars{v}),3) == 12
               LME.(vars{v})(~LME.idxspc(:,:,1:12)) = NaN;
           else
               LME.(vars{v})(~LME.idxspc(:,:,1)) = NaN;
           end
       end
    end

    % calculate area of LME
    LME.area_tot_km2 = ...
        sum(LME.area_km2(LME.idxspc(:,:,1)).*...
        LME.percent_sea(LME.idxspc(:,:,1)),...
        'omitnan');
    disp(['Area = ' num2str(round(LME.area_tot_km2),0) 'km^2']);

    % clean up
    clear tmp_lon vars v

    % Detrend gridded pCO2 using domain mean:
    % Define area weights
    area_weights = LME.area_km2;
    area_weights = repmat(area_weights,1,1,size(LME.fco2_ave_wtd,3));
    area_weights(isnan(LME.fco2_ave_wtd)) = NaN;
    % Calculate area-weighted domain mean
    LME.fco2_dom_mean = ...
        squeeze(sum(sum(LME.fco2_ave_wtd.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    LME.fco2_dom_mean(LME.fco2_dom_mean == 0) = NaN;
    % Fit trend to area weighted domain mean
    yf = leastsq2(LME.month,...
        LME.fco2_dom_mean,0,0,0);
    % Remove difference from mean for each month
    for m = 1:length(LME.month)
        LME.fco2_ave_wtd_detrend(:,:,m) = ...
            LME.fco2_ave_wtd(:,:,m) + ...
            (mean(yf,'omitnan') - yf(m,:));
    end
    % Calculate area-weighted detrended domain mean
    LME.fco2_dom_mean_detrend = ...
        squeeze(sum(sum(LME.fco2_ave_wtd_detrend.*...
        area_weights,1,'omitnan'),2,'omitnan'))./...
        squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
    LME.fco2_dom_mean_detrend(LME.fco2_dom_mean_detrend == 0) = NaN;

    % save gridded pco2 and predictor data in individual LMEs
    if ~isfolder('Data/LME_Data'); mkdir('Data/LME_Data'); end
    save(['Data/LME_Data/' vrs '_' region{n}],'LME','-v7.3');
    clear LME

end

end
