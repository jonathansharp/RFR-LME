% 'CCS' 'GoAK' 'Aleut' 'EBer' 'Beauf' 'BerChuk' 'NEast' 'SEast' 'GoM' 'Car'
% 'Hawaii' 'AmSamoa' 'Jarvis' 'LineIs' 'HowBak' 'Johnst' 'Wake' 'Guam'

options = statset('MaxIter',1000,'Display','final');
Sigma = {'full'};
SharedCovariance = {true};
% num_groups = [5:7;9:11;1:3;1:3;1:3;3:5;1:3;5:7;2:4;2:4;1:3;...
%       ones(1,3);ones(1,3);3:5;ones(1,3);ones(1,3);ones(1,3);ones(1,3)];
num_groups = [6;10;2;2;2;4;2;6;3;3;2;1;1;4;1;1;1;1];
RegularizationValue = 0;
test_idx = 0;