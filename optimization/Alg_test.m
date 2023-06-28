% wrapper to test RFR options 

% Set options
nTrees_vec = 100;
num_nTrees = length(nTrees_vec);
minLeafSize_vec = [1 2 3 5 7 10 15 20];
num_minLeafSize = length(minLeafSize_vec);
numpredictors_vec = 7; % can be any number, this is defined intrnsically
num_numpredictors = length(numpredictors_vec);

% define regions for loop
define_regions_eiwg

% pre-allocate
RMSE = nan(num_nTrees,num_minLeafSize,num_numpredictors);

% for each region
for n = 1:length(region)

% cluster and fit algs
for i = 1:num_nTrees
    for j = 1:num_minLeafSize
        for k = 1:num_numpredictors

            %% define and select region
            define_regions_eiwg; % define regions
            region = region(n); % select region for test step

            %% Set RFR options
            nTrees = nTrees_vec(i);
            minLeafSize = minLeafSize_vec(j);
            numpredictors = numpredictors_vec(k);

            %% Set GMM options and test indices
            set_gmm_options
            num_groups = num_groups(n); % select number of groups for test
            rfr_test_idx = 1; % set RFR index to test

            %% Fit test algorithms
            fit_algs_probs;

            %% Log results
            load(['Data/' region{1} '/us_lme_model_evals_test'],'Val');
            RMSE(i,j,k) = Val.(region{1}).rmse_rfr(end);
            clear Val

        end
    end
end

% plot RMSE
figure; hold on;
bar(1:length(minLeafSize_vec),reshape(RMSE,...
    num_minLeafSize*num_nTrees,num_numpredictors)');
ylim([min(RMSE(:))-1 max(RMSE(:))+1]);
title('RMSE for various MinLeafSize');
xticks(1:length(minLeafSize_vec));
xticklabels({'1' '2' '3' '5' '7' '10' '15' '20'});
xlabel('Minimum Leaf Size');
ylabel('RMSE');
exportgraphics(gcf,['Figures/' region{n} '/RMSE_test_algs.png']);
close

end