% 1. CCS
% 2. GA
% 3. AI
% 4. EBS
% 5. BS
% 6. NBCS
% 7. NE
% 8. SE
% 9. GM
% 10. CS
% 11. PI

options = statset('MaxIter',1000,'Display','off');
Sigma = 'full'; % predictor variabes are correlated (a priori)
SharedCovariance = false; % each cluster has unique covariance shape
num_groups = [5;3;4;4;2;2;2;5;3;2;5]; %  regional clusters by examining RMSE, BIC, and silhouette scores
RegularizationValue = 0;
test_idx = 0;