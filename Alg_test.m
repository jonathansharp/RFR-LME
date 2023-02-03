% wrapper to test RFR options 

% Set options
nTrees_vec = 100;
num_nTrees = length(nTrees_vec);
minLeafSize_vec = [1 2 3 5 8 12 20];
num_minLeafSize = length(minLeafSize_vec);
numpredictors_vec = 2:11;
num_numpredictors = length(numpredictors_vec);

% pre-allocate
RMSE = nan(num_nTrees,num_minLeafSize,num_numpredictors);

% cluster and fit algs
for i = 1:num_nTrees
    for j = 1:num_minLeafSize
        for k = 1:num_numpredictors
            nTrees = nTrees_vec(i);
            minLeafSize = minLeafSize_vec(j);
            numpredictors = numpredictors_vec(k);
            fit_algs;
            load(['Data/' region{1} '/us_lme_models'],'Mods');
            RMSE(i,j,k) = Mods.(region{1}).rmse_rfr(end);
            clear Mods
        end
    end
end

% plot RMSE
figure; hold on;
bar(numpredictors_vec,reshape(RMSE,...
    num_minLeafSize*num_nTrees,num_numpredictors)');
ylim([min(RMSE(:))-1 max(RMSE(:))+1]);
title('RMSE for various MinLeafSize and NumPredictors');
legend({'1' '2' '3' '5' '8' '12' '20'});
xlabel('Number of Predictors');
ylabel('RMSE');
exportgraphics(gcf,['Figures/' region{n} '/RMSE_test_algs.png']);

% find best combination
min(min(RMSE))