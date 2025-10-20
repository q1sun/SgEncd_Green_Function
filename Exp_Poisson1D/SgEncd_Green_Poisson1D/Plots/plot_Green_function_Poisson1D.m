function plot_Green_function_Poisson1D(Green_exact, Green_NN, V_mesh)

format short e

% refined mesh by interpolation
dx = 0:0.002:1;
dy = 0:0.002:1;
[qx,qy] = meshgrid(dx,dy);

%---------------------%
h=figure('NumberTitle','off','Name','Exact GreenFct','Renderer', 'painters', 'Position', [0 0 500 420]);

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',Green_exact');
qz = Ft(qx,qy);
imagesc(qz)

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
cb.Limits = [0 0.25];
set(cb,'TickLabelInterpreter','latex')
axis xy
axis off
set(gca,'FontSize',18);
set(gcf,'color','w');
hold on

x_post = xlabel('$x$','interpreter','latex');
y_post = ylabel('$y$','interpreter','latex');
x_post.Position(2) = h.Position(2) - 10; 
axis on
xticks([1 500])
xticklabels({'0','1'})
yticks([1 500])
yticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-GreenFct-exact.pdf')
%---------------------%


%---------------------%
h=figure('NumberTitle','off','Name','Deep GreenFct','Renderer', 'painters', 'Position', [0 0 500 420]);

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',Green_NN');
qz = Ft(qx,qy);
imagesc(qz)

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
cb.Limits = [0 0.25];
set(cb,'TickLabelInterpreter','latex')
axis xy
axis off
set(gca,'FontSize',18);
set(gcf,'color','w')
hold on

x_post = xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
x_post.Position(2) = h.Position(2) - 10; 
axis on
xticks([1 500])
xticklabels({'0','1'})
yticks([1 500])
yticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-GreenFct-NN.pdf')
%---------------------%


%---------------------%
h=figure('NumberTitle','off','Name','PtErr GreenFct','Renderer', 'painters', 'Position', [0 0 500 420]);

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',Green_exact' - Green_NN');
qz = Ft(qx,qy);
imagesc(qz)

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
set(cb,'TickLabelInterpreter','latex')
axis xy
axis off
set(gca,'FontSize',18);
set(gcf,'color','w')
hold on

x_post = xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
x_post.Position(2) = h.Position(2) - 10; 
axis on
xticks([1 500])
xticklabels({'0','1'})
yticks([1 500])
yticklabels({'0','1'})
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-GreenFct-PtErr.pdf')
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);
length = sqrt(size(V_mesh(1,:), 2));
x = reshape(V_mesh(1,:), length, length);
y = reshape(V_mesh(2,:), length, length);
z = abs(x - y);
color = reshape(Green_NN, size(x));
surf(x, y, z, color);
shading interp
alpha 1

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'FontSize',20);
set(cb,'TickLabelInterpreter','latex')

xlabel('$x$','interpreter','latex','rotation', 45);
ylabel('$y$','interpreter','latex','rotation', -45);
zlabel('$\varphi(x,y)$','interpreter','latex');

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson1D-Green-NN-xyz.pdf')
%---------------------%


end

