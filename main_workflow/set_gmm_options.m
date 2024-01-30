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
%  regional clusters by examining RMSE, BIC, and silhouette scores
% Oct 6th (with ETOPO2 for Bathymetry)
% num_groups = [5;3;4;4;2;2;2;5;3;2;5];
% Dec. 28th (not sure)
% num_groups = [5;4;6;4;4;7;4;5;5;4;3]; 
% Jan 2nd (with ETOPOv2022 and lat/lon to cluster)
% num_groups = [5;5;6;7;4;4;5;4;2;2;5];
% Jan. 11th (with ETOPOv2022 for Bathymetry)
num_groups = [3;4;3;4;3;6;2;5;2;3;3];
RegularizationValue = 0;
gmm_test_idx = 0;
