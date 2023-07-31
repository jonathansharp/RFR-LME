% Log error statistics (no moorings)
% 
% This script loads error statics determined during fCO2 algorithm fits for
% eighteen US large Marine Ecosystems and saves them in a table
% 
% Written by J.D. Sharp: 7/26/23
% Last updated by J.D. Sharp: 7/26/23
% 

%% pre-allocate table
error_stats = table('Size',[length(region) 5],...
    'VariableTypes',{'double' 'double' 'double' 'double' 'double'},...
    'VariableNames',{'Avg. Err.' 'Med. Err.' 'RMSE' 'MAE' 'R2'},...
    'RowNames',region');

%% loop through regions
for n = 1:length(region)

    %% load error statistics
    for en = 1:size(num_groups,2)
        load(['Data/' region{n} '/us_lme_model_evals_no_moorings'],'Val');
        Val_tmp.avg_err_rfr(en) = Val.(region{n}).avg_err_rfr(end);
        Val_tmp.med_err_rfr(en) = Val.(region{n}).med_err_rfr(end);
        Val_tmp.rmse_rfr(en) = Val.(region{n}).rmse_rfr(end);
        Val_tmp.mae_rfr(en) = Val.(region{n}).mae_rfr(end);
        Val_tmp.r2_rfr(en) = Val.(region{n}).r2_rfr(end);
        clear Val
    end

    %% log in table
    error_stats(n,1) = {mean(Val_tmp.avg_err_rfr,'omitnan')};
    error_stats(n,2) = {mean(Val_tmp.med_err_rfr,'omitnan')};
    error_stats(n,3) = {mean(Val_tmp.rmse_rfr,'omitnan')};
    error_stats(n,4) = {mean(Val_tmp.mae_rfr,'omitnan')};
    error_stats(n,5) = {mean(Val_tmp.r2_rfr,'omitnan')};

    %% clean up
    clear Mods

end

%% save table of error statistics
writetable(error_stats,['Data/ErrorStatistics-no-moorings-' date '.xls']);
