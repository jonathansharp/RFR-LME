% 'CCS' 'GoAK' 'Aleut' 'EBer' 'Beauf' 'BerChuk' 'NEast' 'SEast' 'GoM' 'Car'
% 'Hawaii' 'AmSamoa' 'Jarvis' 'LineIs' 'HowBak' 'Johnst' 'Wake' 'Guam'

options = statset('MaxIter',1000,'Display','off');
Sigma = 'full';
SharedCovariance = false;
num_groups = [6;10;2;2;2;4;2;6;3;3;2;1;1;4;1;1;1;1];
RegularizationValue = 0;
test_idx = 0;