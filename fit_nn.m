function [nn,err,rmse,r2,Y_fit,delta] = ...
    fit_nn(X_mod,Y_mod,idx_vars,idx_group)

% set training function
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation

% set hidden layer sizes
size1 = [20 15 10];
size2 = [10 15 20];

% pre-allocate error stats
err = nan(max(idx_group),length(size1));
rmse = nan(max(idx_group),length(size1));
r2 = nan(max(idx_group),length(size1));

% for each cluster
for c = 1:max(idx_group)

    % cluster label
    clab = ['c' num2str(c)];
    
    % define X and Y datasets for each cluster
    X_temp = X_mod(idx_group==c,idx_vars);
    Y_temp = Y_mod(idx_group==c);

    % pre-allocate test vectors
    Y_fit.(clab) = nan(sz,length(size1));
    delta.(clab) = nan(sz,length(size1));

    % for each architecture
    for a = 1:length(size1)

        % architecture label
        alab = ['a' num2str(a)];
    
        % create net
        net.(clab).(alab) = feedforwardnet([size1(a) size2(a)],trainFcn);
    
        % set up division of data for training, validation, testing
        net.(clab).(alab).divideParam.trainRatio = 85/100; % default: 70/100
        net.(clab).(alab).divideParam.valRatio = 15/100; % default: 15/100
        net.(clab).(alab).divideParam.testRatio = 0/100;% default: 15/100
        
        % set up training parameters
        net.(clab).(alab).trainParam.max_fail = 6; % default: 6
        net.(clab).(alab).trainParam.mu_max = 1e10; % default: 1e10
        net.(clab).(alab).trainParam.min_grad = 1e-7; % default: 1e-7
    
        % set up cross-validation
        numFolds = 5;
        sz = length(Y_temp);
        ints = randperm(length(Y_temp))';

        % for each fold
        for f = 1:numFolds

            % fold label
            flab = ['f' num2str(f)];
        
            % index for cross-validation training and test
            idx_test = ints > (f-1) * (sz/numFolds) & ints <= f * (sz/numFolds);
            idx_train = ~idx_test;
            
            % get train and test data for this fold
            X_train = X_temp(idx_train,:);
            Y_train = Y_temp(idx_train);
            X_test = X_temp(idx_test,:);
            Y_test = Y_temp(idx_test);

            % train network
            nn_cv.(clab).(alab).(flab) = ...
                train(net.(clab).(alab),X_train',Y_train');
            
            % test network
            net_temp = nn_cv.(clab).(alab).(flab);
            Y_fit.(clab)(idx_test,a) = net_temp(X_test')';
            delta.(clab)(idx_test,a) = Y_fit.(clab)(idx_test,a)-Y_test;
        
        end

        % save network error statistics for each fold
        err(c,a) = mean(delta.(clab)(:,a));
        rmse(c,a) = sqrt(mean(delta.(clab)(:,a).^2));
        r2(c,a) = corr(Y_fit.(clab)(:,a),Y_temp).^2;

        % train network with all data
        nn.(clab).(alab) = train(net.(clab).(alab),X_temp',Y_temp');

    end

    % average each architecture
    Y_fit.(clab)(:,a+1) = mean(Y_fit.(clab),2);
    delta.(clab)(:,a+1) = mean(delta.(clab),2);

    % save network error statistics
    err(c,a+1) = mean(delta.(clab)(:,a+1));
    rmse(c,a+1) = sqrt(mean(delta.(clab)(:,a+1).^2));
    r2(c,a+1) = corr(Y_fit.(clab)(:,a+1),Y_temp).^2;

end

% end function
end
