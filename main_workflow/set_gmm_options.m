% 'CCS' 'GoAK' 'Aleut' 'EBer' 'Beauf' 'BerChuk' 'NEast' 'SEast' 'GoM' 'Car'
% 'Hawaii' 'AmSamoa' 'Jarvis' 'LineIs' 'HowBak' 'Johnst' 'Wake' 'Guam'

options = statset('MaxIter',1000,'Display','off');
Sigma = 'full'; % predictor variabes are correlated (a priori)
SharedCovariance = false; % each cluster has unique covariance shape
num_groups = [5;7;3;4;4;3;3;3;3;5;2;... % large regional clusters determined by examining RMSE, BIC, and silhouette scores
    1;1;1;1;1;1;1]; % island regions all one cluster
RegularizationValue = 0;
test_idx = 0;