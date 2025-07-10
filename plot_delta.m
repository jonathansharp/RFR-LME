% Plot delta pCO2
% 
% Written by J.D. Sharp: 9/28/22
% Last updated by J.D. Sharp: 11/29/22
% 

function plot_delta(x_edges,x_mids,y_edges,y_mids,Y_fit,delta,reg,xlimit,ylimit,type)

% assemble histogram counts
counts = histcounts2(Y_fit,delta,x_edges,y_edges);

% plot figure
figure; hold on; box off;
set(gcf,'Position',[100 100 800 400]);
h=imagesc(x_mids,y_mids,counts');
plot([min(x_edges) max(x_edges)],[0 0],'k--');
xlim(xlimit); ylim(ylimit);
xlabel('{\itp}CO_{2} (\muatm)');
ylabel('\Delta{\itp}CO_{2} (\muatm)');
text(max(xlimit)-200,min(ylimit)+50,['RMSE = ' ...
    num2str(round(sqrt(mean(delta.^2,'omitnan')),2))],'fontsize',16);
myColorMap = parula(20);
myColorMap(1,:) = 1;
colormap(gca,myColorMap);
set(gca,'ColorScale','log');
max_c = numel(num2str(max(counts(:))));
caxis(gca,[1e0 1*10^max_c]);
c=colorbar;
c.TickLength = 0;
c.Label.String = 'Frequency';
c.Label.FontSize = 14;

% save figure
if ~isfolder('Figures'); mkdir('Figures'); end
exportgraphics(gcf,['Figures/' reg '_del_fCO2_' type '.png']);
close
