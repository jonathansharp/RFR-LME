function log_errs(vrs,num_groups,region)

%% pre-allocate table
error_stats = table('Size',[length(region) 7],...
    'VariableTypes',{'double' 'double' 'double' 'double' 'double' 'double' 'double'},...
    'VariableNames',{'Avg. Err.' 'Avg. Abs. Err.' 'Med. Err.' 'Med. Abs. Err' 'IQR' 'RMSE' 'R2'},...
    'RowNames',region);

%% loop through regions
for n = 1:length(region)

    %% load error statistics and calculate some extra ones
    load(['Models/' region{n} '/us_lme_models_' vrs '_c' num2str(num_groups(n))]);
    rfr.avg_abs_err = mean(abs(rfr.delta.all(:,end)),'omitnan');
    rfr.med_abs_err = median(abs(rfr.delta.all(:,end)),'omitnan');
    rfr.iqr = iqr(rfr.delta.all(:,end));

    %% log in table
    error_stats(n,1) = {rfr.avg_err(end)};
    error_stats(n,2) = {mean(rfr.avg_abs_err,'omitnan')};
    error_stats(n,3) = {rfr.med_err(end)};
    error_stats(n,4) = {mean(rfr.med_abs_err,'omitnan')};
    error_stats(n,5) = {mean(rfr.iqr,'omitnan')};
    error_stats(n,6) = {rfr.rmse(end)};
    error_stats(n,7) = {rfr.r2(end)};

end

%% save table of error statistics
writetable(error_stats,['IndsAndStats/ErrorStatistics-' date '.xls']);
