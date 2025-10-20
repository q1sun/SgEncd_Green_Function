function plot_additive_prectd_eigvalues(lambda_max, lambda_min, index)

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot([lambda_min,lambda_max],0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0, 3]);
ylim([0,0.1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,lambda_min, 1, 2, lambda_max, 3]);
set(gca,'XTickLabel',{'$0$','$\lambda_{\mathrm{min}}$','$1$','2', '$\lambda_{\mathrm{max}}$','$3$'});
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, sprintf('Figure/Meshsize_h/fig-Poisson1D-eigenval-green-h%d.pdf',index));
%---------------------%

end