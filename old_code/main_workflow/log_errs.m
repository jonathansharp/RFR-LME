% Log error statistics
% 
% This script loads error statics determined during fCO2 algorithm fits for
% eighteen US large Marine Ecosystems and saves them in a table
% 
% Written by J.D. Sharp: 1/10/23
% Last updated by J.D. Sharp: 1/10/23
% 

%% pre-allocate table
error_stats = table('Size',[length(region) 7],...
    'VariableTypes',{'double' 'double' 'double' 'double' 'double' 'double' 'double'},...
    'VariableNames',{'Avg. Err.' 'Avg. Abs. Err.' 'Med. Err.' 'Med. Abs. Err' 'IQR' 'RMSE' 'R2'},...
    'RowNames',region);

%% loop through regions
for n = 1:length(region)

    %% load error statistics
    for en = 1:size(num_groups,2)
        load(['Data/' region{n} '/us_lme_model_evals'],'Val');
        Val_tmp.avg_err_rfr(en) = mean(Val.(region{n}).delta_rfr.all(:,end));
        Val_tmp.avg_abs_err_rfr(en) = mean(abs(Val.(region{n}).delta_rfr.all(:,end)));
        Val_tmp.med_err_rfr(en) = median(Val.(region{n}).delta_rfr.all(:,end));
        Val_tmp.med_abs_err_rfr(en) = median(abs(Val.(region{n}).delta_rfr.all(:,end)));
        Val_tmp.iqr_rfr(en) = iqr(Val.(region{n}).delta_rfr.all(:,end));
        Val_tmp.rmse_rfr(en) = sqrt(mean(Val.(region{n}).delta_rfr.all(:,end).^2));
        Val_tmp.r2_rfr(en) = Val.(region{n}).r2_rfr(end);
        clear Val
    end

    %% log in table
    error_stats(n,1) = {mean(Val_tmp.avg_err_rfr,'omitnan')};
    error_stats(n,2) = {mean(Val_tmp.avg_abs_err_rfr,'omitnan')};
    error_stats(n,3) = {mean(Val_tmp.med_err_rfr,'omitnan')};
    error_stats(n,4) = {mean(Val_tmp.med_abs_err_rfr,'omitnan')};
    error_stats(n,5) = {mean(Val_tmp.iqr_rfr,'omitnan')};
    error_stats(n,6) = {mean(Val_tmp.rmse_rfr,'omitnan')};
    error_stats(n,7) = {mean(Val_tmp.r2_rfr,'omitnan')};

    %% clean up
    clear Mods

end

%% save table of error statistics
writetable(error_stats,['IndsAndStats/ErrorStatistics-' date '.xls']);
