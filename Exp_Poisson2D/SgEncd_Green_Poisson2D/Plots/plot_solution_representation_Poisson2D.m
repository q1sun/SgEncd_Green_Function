function plot_solution_representation_Poisson2D(V_mesh, u_exact_MeshPts, u_repNG_MeshPts)

format short e

dx = -1:.002:1;
dy = -1:.002:1;
[qx,qy] = meshgrid(dx,dy);

for i = 1 : size(qx,1)
    for j = 1 : size(qx,2)
        if sqrt((qx(i,j)).^2+(qy(i,j)).^2) > 1 + 0.0001
            qx(i,j) = nan;
            qy(i,j) = nan;
        end
    end
end


%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','exact solution on unit disc');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',u_exact_MeshPts(:));
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
exportgraphics(ax, 'Figures/Task2/fig-Poisson2D-u-exact.pdf')
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','represented solution on unit disc');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',u_repNG_MeshPts(:));
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
exportgraphics(ax, 'Figures/Task2/fig-Poisson2D-u-repNG.pdf')
%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%
figure('NumberTitle','off','Name','pointwise error for represented solution');

Ft = TriScatteredInterp(V_mesh(1,:)',V_mesh(2,:)',u_exact_MeshPts(:) - u_repNG_MeshPts(:));
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
exportgraphics(ax, 'Figures/Task2/fig-Poisson2D-u-PtErr.pdf')
%%%%%%%%%%%%%%%


end
