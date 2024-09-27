% This function plots GMM test parameters
% 
% Written by J.D. Sharp: 5/2/23
% Last updated by J.D. Sharp: 5/2/23
% 

% define regions for loop
define_regions

% initiate figure
t=tiledlayout(6,3);
t.Padding = 'Compact';
t.TileSpacing = 'Compact';
clrs = colororder;

xlabel('No. Clusters');
    ylabel('Root Mean Squared Error (\mumol kg^{-1})');
    ylabel('Bayesian Information Criterion');

for n = 1:length(region)

    %% Set GMM options
    set_gmm_options % set baseline options
    num_groups_vec = 1:20; % test different group amounts
    nGrp = length(num_groups_vec);
    test_idx = 1; % set index to test

    % load error stats
    load(['Data/' region{n} '/RMSE'],'RMSE');
    load(['Data/' region{n} '/BIC'],'BIC');
    load(['Data/' region{n} '/SIL'],'SIL');

    % plot combined
    ax=nexttile; hold on;
    yyaxis left
    plot(num_groups_vec(2:end),RMSE(2:end),'color',clrs(1,:),'linewidth',2);
    plot([num_groups(n) num_groups(n)],[min(RMSE(:))-1 max(RMSE(:))+1],'k--');
    ylim([min(RMSE(:))-0.5 max(RMSE(:))+0.5]);
    yyaxis right
    plot(num_groups_vec(2:end),BIC(2:end),'color',clrs(2,:),'linewidth',3);
    set(gca,'YDir','reverse','Color','none','XColor','none','YColor','none');
    xlim([1.5 20.5]);
    plot(num_groups_vec(2:end),SIL(2:end),'color',clrs(4,:),'linewidth',3);
    set(gca,'Color','none','YDir','reverse','XColor','none','YColor',clrs(4,:));
    plot(num_groups_vec(2:end),SIL(2:end),'color',clrs(4,:),'linestyle','none');
    ylabel('Global Mean Silhouette Score');
    exportgraphics(gcf,['Figures/' region{1} '_Combined_test_groups.png']);
    close