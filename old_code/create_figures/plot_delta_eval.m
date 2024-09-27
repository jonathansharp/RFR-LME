% Plot 2D histogram of discrete measurements compared to gridded values
% 
% Written by J.D. Sharp: 2/17/23
% Last updated by J.D. Sharp: 2/17/23
% 

function plot_delta_eval(edges,discrete,gridded,var_type,var_lab,units,rounder,source)

% assemble histogram counts
counts = histcounts2(discrete,gridded,edges,edges);

% plot figure
figure; hold on; box off;
mids = diff(edges) + edges(1:end-1);
h=imagesc(mids,mids,counts');
plot([min(edges) max(edges)],[min(edges) max(edges)],'k--');
xlabel([var_type ' (' source ')']);
ylabel([var_type ' (LME-RFR)']);
xlim([min(edges) max(edges)]);
ylim([min(edges) max(edges)]);
myColorMap = parula(40);
myColorMap(1,:) = 1;
colormap(gca,myColorMap);
max_c = round(max(counts(:)),1);
caxis(gca,[0 max_c]);
c=colorbar;
c.TickLength = 0;
c.Label.String = 'Frequency';
c.Label.FontSize = 14;
rmse = sqrt(mean((discrete - gridded).^2,'omitnan'));
text(edges(end-15),edges(5),['RMSE = ' num2str(round(rmse,rounder))]);

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,['Figures/del_' var_type '_' source '.png']);
close
