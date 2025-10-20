function plot_fixedY_Green_function_Poisson2D(exact_G_fixedY1, NN_G_fixedY1, exact_G_fixedY2, NN_G_fixedY2, V_mesh)

format short e

dx = -1:.01:1;
dy = -1:.01:1;
[qx,qy] = meshgrid(dx,dy);

for i = 1 : size(qx,1)
    for j = 1 : size(qx,2)
        if sqrt((qx(i,j)).^2+(qy(i,j)).^2) > 1 + 0.0001
            qx(i,j) = nan;
            qy(i,j) = nan;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',NN_G_fixedY1(:));
qz = Ft(qx,qy);
[nr,nc] = size(qz);

surf(qx,qy, qz);
shading interp
alpha 0.5

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'FontSize',18);
set(cb,'TickLabelInterpreter','latex')

xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$\check{G}(x_1,x_2,0,0)$','interpreter','latex');

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-SurfGreen-NN-y1.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('NumberTitle','off','Name','Higher Dimensional Func','Renderer', 'painters', 'Position', [0 0 500 420]);

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',NN_G_fixedY2(:));
qz = Ft(qx,qy);
[nr,nc] = size(qz);

surf(qx,qy, qz);
shading interp
alpha 0.5

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'FontSize',18);
set(cb,'TickLabelInterpreter','latex')

xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$\check{G}(x_1,x_2,0.7,0.2)$','interpreter','latex');

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-SurfGreen-NN-y2.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','deep Green on unit disc');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',NN_G_fixedY1(:));
qz = Ft(qx,qy);
[nr,nc] = size(qz);
pcolor([qz nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(cb,'TickLabelInterpreter','latex')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'dataAspectRatio',[1 1 1])
set(gca,'FontSize',18);
set(gca,'xtick',[1 length(dx)],'xticklabel',[-1 1])
set(gca,'ytick',[1 length(dy)],'yticklabel',[-1 1])
axis xy
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gcf,'color','w')

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-GreenFct-NN-y1.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','pointwise error Green');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)', abs(exact_G_fixedY1(:) - NN_G_fixedY1(:)));
qz = Ft(qx,qy);
[nr,nc] = size(qz);
pcolor([qz nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(cb,'TickLabelInterpreter','latex')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'dataAspectRatio',[1 1 1])
set(gca,'FontSize',18);
set(gca,'xtick',[1 length(dx)],'xticklabel',[-1 1])
set(gca,'ytick',[1 length(dy)],'yticklabel',[-1 1])
axis xy
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gcf,'color','w')

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-GreenFct-PtErr-y1.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','deep Green on unit disc');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',NN_G_fixedY2(:));
qz = Ft(qx,qy);
[nr,nc] = size(qz);
pcolor([qz nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(cb,'TickLabelInterpreter','latex')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'dataAspectRatio',[1 1 1])
set(gca,'FontSize',18);
set(gca,'xtick',[1 length(dx)],'xticklabel',[-1 1])
set(gca,'ytick',[1 length(dy)],'yticklabel',[-1 1])
axis xy
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gcf,'color','w')

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-GreenFct-NN-y2.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','pointwise error Green');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)', abs(exact_G_fixedY2(:) - NN_G_fixedY2(:)));
qz = Ft(qx,qy);
[nr,nc] = size(qz);
pcolor([qz nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');

colormap jet
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(cb,'TickLabelInterpreter','latex')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'dataAspectRatio',[1 1 1])
set(gca,'FontSize',18);
set(gca,'xtick',[1 length(dx)],'xticklabel',[-1 1])
set(gca,'ytick',[1 length(dy)],'yticklabel',[-1 1])
axis xy
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gcf,'color','w')

ax = gcf;
exportgraphics(ax, 'Figures/Task1/fig-Poisson2D-GreenFct-PtErr-y2.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end