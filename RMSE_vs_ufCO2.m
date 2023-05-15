
% this script defines the bounds of the eighteen LMEs
define_regions

for n = 1:length(region)

    %% load estimated OA grid and model evaluation
    load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    load(['Data/' region{n} '/us_lme_model_evals'],'Val');

    %% print uncertainties
    disp([region{n} ', RMSE = ' num2str(round(Val.(region{n}).rmse_rfr(end),1))]);
    disp([region{n} ', ufCO2 = ' num2str(round(mean(OAI_grid.(region{n}).ufCO2(:),'omitnan'),1))]);

end