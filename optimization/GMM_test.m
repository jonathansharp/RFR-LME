% Wrapper to test GMM 
% 
% This function tests the number of clusters to include in the GMM
% 
% Written by J.D. Sharp: 9/20/22
% Last updated by J.D. Sharp: 1/2/24
% 

addpath(genpath(pwd));

% define regions for loop
define_regions_eiwg

for n = 1:length(region)

    %% define and select region
    define_regions_eiwg; % define regions
    region = region(n); % select region for test step

    %% Set GMM options
    set_gmm_options % set baseline options
    num_groups_vec = 1:20; % test different group amounts
    nGrp = length(num_groups_vec);
    gmm_test_idx = 1; % set GMM index to test

    %% Set RFR options
    set_rfr_options % set baseline options
    nTrees = 100; % reduce number of trees for efficiency
    rfr_test_idx = 1;
    
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
    title(['RMSE for various {\itk} Choices (' Sigma ')']);
    xlabel('Number of Groups');
    ylabel('RMSE');
    text(3,max(BIC(:)),['Best: ' num2str(num_groups_vec(Ri))]);
    exportgraphics(gcf,['Figures/' region{1} '_RMSE_test_groups.png']);
    close

    % plot BIC
    figure; hold on;
    bar(num_groups_vec,BIC);
    title(['BIC for various {\itk} Choices (' Sigma ')']);
    % legend({'diagonal, shared' 'diagonal, unshared' 'full, shared' 'full, unshared'});
    xlabel('Number of Groups ({\itk})');
    ylabel('BIC');
    text(3,max(BIC(:)),['Best: ' num2str(num_groups_vec(Bi))]);
    exportgraphics(gcf,['Figures/' region{1} '_BIC_test_groups.png']);
    close

    % plot mean SIL
    figure; hold on;
    bar(num_groups_vec,SIL);
    title(['SIL for various {\itk} Choices (' Sigma ')']);
    % legend({'diagonal, shared' 'diagonal, unshared' 'full, shared' 'full, unshared'});
    xlabel('Number of Groups ({\itk})');
    ylabel('SIL');
    text(3,max(SIL(:)),['Best: ' num2str(num_groups_vec(Si))]);
    exportgraphics(gcf,['Figures/' region{1} '_SIL_test_groups.png']);
    close

    % plot combined
    clrs = colororder;
    set_gmm_options;
    figure; hold on;
    set(gcf,'Position',[720 478 700 420])
    set(gca,'Position',[0.2300 0.1100 0.6750 0.8150]);
    xlim([1.5 20.5]);
    yyaxis left
    plot(num_groups_vec(2:end),RMSE(2:end),'color',clrs(1,:),'linewidth',3);
    %plot([num_groups(n) num_groups(n)],[min(RMSE(:))-1 max(RMSE(:))+1],'k--');
    ylim([min(RMSE(:))-0.5 max(RMSE(:))+0.5]);
    xlabel('No. Clusters');
    ylabel('Root Mean Squared Error (\mumol kg^{-1})');
    yyaxis right
    plot(num_groups_vec(2:end),BIC(2:end),'color',clrs(2,:),'linewidth',3);
    ylabel('Bayesian Information Criterion');
    axes('position',[0.2300 0.1100 0.6750 0.8150]); hold on
    set(gca,'YDir','normal','Color','none','XColor','none','YColor','none');
    xlim([1.5 20.5]);
    plot(num_groups_vec(2:end),SIL(2:end),'color',clrs(4,:),'linewidth',3);
    axes('Position',[0.1100 0.1100 0.6750 0.8150]); hold on
    set(gca,'Color','none','YDir','normal','XColor','none','YColor',clrs(4,:));
    plot(num_groups_vec(2:end),SIL(2:end),'color',clrs(4,:),'linestyle','none');
    ylabel('Global Mean Silhouette Score');
    exportgraphics(gcf,['Figures/' region{1} '_Combined_test_groups.png']);
    close

    % save error stats
    save(['Data/' region{1} '/RMSE'],'RMSE');
    save(['Data/' region{1} '/BIC'],'BIC');
    save(['Data/' region{1} '/SIL'],'SIL');

    % clean up
    clear RMSE BIC region

end

