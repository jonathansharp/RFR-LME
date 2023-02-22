% Wrapper to test GMM 
% 
% This function tests the number of clusters to include in the GMM
% 
% Written by J.D. Sharp: 9/20/22
% Last updated by J.D. Sharp: 2/21/23
% 

% define regions for loop
define_regions

for n = 1:length(region)

    %% define and select region
    define_regions; % define regions
    region = region(n); % select region for test step

    %% Set GMM options
    set_gmm_options % set baseline options
    num_groups_vec = [1:20]; % test different group amounts
    nGrp = length(num_groups_vec);
    test_idx = 1; % set index to test

    %% Set RFR options
    set_rfr_options % set baseline options
    nTrees = 100; % reduce number of trees
    
    %% pre-allocate
    RMSE = nan(nGrp,1);
    BIC = nan(nGrp,1);
    SIL = nan(nGrp,1);
    
    %% cluster and fit algs
    for k = 1:nGrp
        num_groups = num_groups_vec(k);
        load(['Data/' region{1} '/variable_arrays'],'Vars_array');
        % Only attempt test for cluster amounts that are
        % appropriate for the amount of training data
        if num_groups < length(Vars_array.(region{1}).Y_mod)/100
            cluster_on_grid; % cluster
            fit_algs_probs; % fit algorithm
            load(['Data/' region{1} '/us_lme_model_evals_test'],'Val');
            RMSE(k) = Val.(region{1}).rmse_rfr(end);
            load(['Data/' region{1} '/gridded_clusters_test'],'Clusts_grid');
            BIC(k) = Clusts_grid.(region{1}).BIC;
            SIL(k) = mean(Clusts_grid.(region{1}).SIL);
            clear Val Clusts_grid
        else
            RMSE(k) = NaN;
            BIC(k) = NaN;
            SIL(k) = NaN;
        end
    end
    
    %% find number of groups for minimum SIL

     Ri = find(RMSE==min(RMSE(:)));
     Bi = find(BIC==min(BIC(:)));
     Si = find(SIL==max(SIL(:)));


    % plot RMSE
    figure; hold on;
    bar(num_groups_vec,RMSE);
    ylim([min(RMSE(:))-1 max(RMSE(:))+1]);
    title('RMSE for various {\itk} Choices');
    xlabel('Number of Groups');
    ylabel('RMSE');
    text(3,max(BIC(:)),['Best: ' num2str(num_groups_vec(Ri))]);
    exportgraphics(gcf,['Figures/' region{1} '_RMSE_test_groups.png']);
    close

    % plot BIC
    figure; hold on;
    bar(num_groups_vec,BIC);
    title('BIC for various {\itk} Choices');
    % legend({'diagonal, shared' 'diagonal, unshared' 'full, shared' 'full, unshared'});
    xlabel('Number of Groups ({\itk})');
    ylabel('BIC');
    text(3,max(BIC(:)),['Best: ' num2str(num_groups_vec(Bi))]);
    exportgraphics(gcf,['Figures/' region{1} '_BIC_test_groups.png']);
    close

    % plot mean SIL
    figure; hold on;
    bar(num_groups_vec,SIL);
    title('SIL for various {\itk} Choices');
    % legend({'diagonal, shared' 'diagonal, unshared' 'full, shared' 'full, unshared'});
    xlabel('Number of Groups ({\itk})');
    ylabel('SIL');
    text(3,max(SIL(:)),['Best: ' num2str(num_groups_vec(Si))]);
    exportgraphics(gcf,['Figures/' region{1} '_SIL_test_groups.png']);
    close

    % save error stats
    save(['Data/' region{1} '/RMSE'],'RMSE');
    save(['Data/' region{1} '/BIC'],'BIC');

    % clean up
    clear RMSE BIC region

end