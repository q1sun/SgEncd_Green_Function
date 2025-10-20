function plot_solution_representation_Poisson1D( V_mesh, uh_exact, uh_NN, fh_exact )

format short e

%---------------------%
h=figure('NumberTitle','off','Name','Exact Solution','Renderer', 'painters', 'Position', [0 0 500 420]);
plot( V_mesh(:), fh_exact(:), 'b-','linewidth',1)

xlabel('$x$','interpreter','latex');
ylabel('forcing term','interpreter','latex');
legend('$f(x)$','interpreter','latex','Location','south');
axis xy
set(gca,'FontSize',20);
set(gcf,'color','w')
xticks([0 1])
xticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task2/fig-Poisson1D-Source-f.pdf')
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Exact Solution','Renderer', 'painters', 'Position', [0 0 500 420]);
plot( V_mesh(:), uh_NN(:), 'r--','linewidth',2)
hold on
plot( V_mesh(:), uh_exact(:), 'k-','linewidth',1)

xlabel('$x$','interpreter','latex');
ylabel('solution value','interpreter','latex');
legend('$\check{u}(x)$','$u(x)$','interpreter','latex','Location','south');
axis xy
set(gca,'FontSize',20);
set(gcf,'color','w')
xticks([0 1])
xticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task2/fig-Poisson1D-Solution.pdf')
%---------------------%



%---------------------%
h=figure('NumberTitle','off','Name','PtErr for Solution Prediction','Renderer', 'painters', 'Position', [0 0 500 420]);
plot( V_mesh(:), log(abs(uh_exact(:) - uh_NN(:))), 'color', [0.4660 0.6740 0.1880] ,'linewidth',1.2)

xlabel('$x$','interpreter','latex');
ylabel('error in log10-scale','interpreter','latex');
legend('log$_{10}|u(x)-\check{u}(x)|$','interpreter','latex','Location','south');
axis xy
set(gca,'FontSize',20);
set(gcf,'color','w')
xticks([0 1])
xticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

bx = gca;
% bx.YAxis.Exponent = -3;
ax = gcf;
% saveas(ax, 'Figures/Task2/fig-Poisson1D-Solution-PtErr.pdf');
exportgraphics(ax, 'Figures/Task2/fig-Poisson1D-Solution-PtErr.pdf')
%---------------------%


