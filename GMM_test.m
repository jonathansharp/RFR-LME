% Wrapper to test GMM 
% 
% 
% 
% Written by J.D. Sharp: 9/20/22
% Last updated by J.D. Sharp: 1/12/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions

for n = 1:length(region)

    % define and select region
    define_regions;
    region = region(n);

    % Set GMM options
    options = statset('MaxIter',1000);
    Sigma_vec = {'diagonal' 'full'};
    nSigma = numel(Sigma_vec);
    SharedCovariance_vec = {true false};
    SCtext = {'shared' 'unshared'};
    nSC = numel(SharedCovariance_vec);
    num_groups_vec = [1:20];
    nGrp = length(num_groups_vec);
    RegularizationValue = 0;
    test_idx = 1;
    
    % pre-allocate
    RMSE = nan(nSigma,nSC,nGrp);
    BIC = nan(nSigma,nSC,nGrp);
    
    % cluster and fit algs
    for i = 1:nSigma
        for j = 1:nSC
            for k = 1:nGrp
                Sigma = Sigma_vec(i);
                SharedCovariance = SharedCovariance_vec(j);
                num_groups = num_groups_vec(k);
                % Only attempt test for cluster amounts that are
                % appropriate for the amount of training data
                load(['Data/' region{1} '/variable_arrays'],'Vars_array');
                if num_groups < length(Vars_array.(region{1}).Y_mod)/100
                    cluster_on_grid;
                    set_rfr_options;
                    fit_algs;
                    load(['Data/' region{1} '/us_lme_model_evals_test'],'Val');
                    RMSE(i,j,k) = Val.(region{1}).rmse_rfr(end);
                    load(['Data/' region{1} '/gridded_clusters_test'],'Clusts_grid');
                    BIC(i,j,k) = Clusts_grid.(region{1}).BIC;
                    clear Val Clusts_grid
                else
                    RMSE(i,j,k) = NaN;
                    BIC(i,j,k) = NaN;
                end
            end
        end
    end
    
    % find configuration for minimum RMSE/BIC
    [iR,jR,kR] = ind2sub(size(RMSE),find(RMSE==min(RMSE(:))));
    [iB,jB,kB] = ind2sub(size(BIC),find(BIC==min(BIC(:))));

    % plot RMSE
    figure; hold on;
    bar(num_groups_vec,reshape(RMSE,nSigma*nSC,nGrp));
    ylim([min(RMSE(:))-1 max(RMSE(:))+1]);
    title('RMSE for various {\itk} and {\it\Sigma} Choices');
    legend({'diagonal, shared' 'full, shared' 'diagonal, unshared' 'full, unshared'});
    xlabel('Number of Groups');
    ylabel('RMSE');
    text(1,max(RMSE(:))+0.8,['Best: ' Sigma_vec{iR} ', ' SCtext{jR} ', ' num2str(num_groups_vec(kR))]);
    exportgraphics(gcf,['Figures/' region{1} '_RMSE_test_groups.png']);
    close

    % plot BIC
    figure; hold on;
    bar(num_groups_vec,reshape(BIC,nSigma*nSC,nGrp)');
    title('BIC for various {\itk} and {\it\Sigma} Choices');
    legend({'diagonal, shared' 'diagonal, unshared' 'full, shared' 'full, unshared'});
    xlabel('Number of Groups ({\itk})');
    ylabel('BIC');
    text(3,max(BIC(:)),['Best: ' Sigma_vec{iB} ', ' SCtext{jB} ', ' num2str(num_groups_vec(kB))]);
    exportgraphics(gcf,['Figures/' region{1} '_BIC_test_groups.png']);
    close

    % save error stats
    save(['Data/' region{1} '/RMSE'],'RMSE');
    save(['Data/' region{1} '/BIC'],'BIC');

    % clean up
    clear RMSE BIC region

end