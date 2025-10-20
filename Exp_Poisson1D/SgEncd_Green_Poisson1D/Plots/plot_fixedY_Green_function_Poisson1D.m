function plot_fixedY_Green_function_Poisson1D(Green_exact_trace, Green_NN, Green_NN_trace, V_mesh, mesh_x)

format short e

% fixed y = 0.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);
length = sqrt(size(V_mesh(1,:),2));
x = reshape(V_mesh(1,:), length, length);
z = reshape(V_mesh(2,:), length, length);
Green_NN = reshape(Green_NN, length, length);

surf(x, z, Green_NN);
shading interp
alpha 0.5

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'FontSize',20);
set(cb,'TickLabelInterpreter','latex')

hold on
x = mesh_x;
z = abs(x - 0.3);
c1 = plot3(x, z, Green_NN_trace,'k','linewidth',3);

view(-30,-10);

xlabel('$x$','interpreter','latex');
ylabel('$\varphi$','interpreter','latex');
zlabel('$\widehat{G}(x,y,z)$','interpreter','latex');
legend(c1,'$\widehat{G}(x,0.3,\varphi(x,0.3))$','interpreter','latex','Position',[0.2,0.455,0.2,1])

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-HigherFunc-y-0.3.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);

set(gca,'FontSize',20);

plot(x, Green_NN_trace,'r--','linewidth',2);
hold on
plot(x, Green_exact_trace, 'k','linewidth',1);

ylim([0 0.25]);

xlabel('$x$','interpreter','latex');
ylabel('function value', 'interpreter','latex');
legend('$\check{G}(x,0.3)$','$G(x,0.3)$','interpreter','latex','Position',[0.56,0.33,0.2,1]);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(gca,'FontSize',20);
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-1DGreenFunc-y-0.3.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);

plot(x,log(abs(Green_exact_trace - Green_NN_trace)),'color', [0.4660 0.6740 0.1880],'linewidth',1.2);

xlabel('$x$','interpreter','latex');
ylabel('error in log10-scale', 'interpreter','latex');
legend('log$_{10}|\check{G}(x,0.3) - G(x,0.3) |$','interpreter','latex','Location','north');

set(gca,'FontSize',20);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-1DGreenFunc-y2-error.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end