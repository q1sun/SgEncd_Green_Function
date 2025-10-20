function plot_FreErr_Helmholtz1D(FreErr_Hybrid,FreErr_Jacobi,index_num,mtd)

iters=size(FreErr_Hybrid,1)-1;

min_E=min(min(FreErr_Hybrid));
max_E=max(max(FreErr_Hybrid));

%---------------------%
h=figure('NumberTitle','off','Name','Magnitude of ModeWise Error','Renderer', 'painters', 'Position', [0 0 700 350]);

iter_index=repmat((1:1:iters+1)',size(FreErr_Jacobi,2),1);
mode_index=reshape(repmat(1:1:size(FreErr_Jacobi,2),size(FreErr_Jacobi,1),1),size(FreErr_Jacobi,1)*size(FreErr_Jacobi,2),1);
Ft = TriScatteredInterp(iter_index,mode_index,reshape(FreErr_Jacobi,size(FreErr_Jacobi,1)*size(FreErr_Jacobi,2),1));
[qx,qy] = meshgrid(1:0.1:iters+1,1:0.1:size(FreErr_Jacobi,2));
qz = Ft(qx,qy);
x=[1 size(FreErr_Jacobi,1)];
y=[1 size(FreErr_Jacobi,2)];

imagesc(x,y,qz)

colormap winter
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'ColorScale','log')
set(cb,'TickLabelInterpreter','latex')
axis xy
axis off
set(gca,'FontSize',18);
set(gcf,'color','w');
hold on
caxis([min_E max_E])

axis on
xticks(1:10:51)
xticklabels({'1','11','21','31','41','51'})
xlabel('index $k$ of iterations','interpreter','latex');
ylabel('index $j$ of frequencies','interpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ax = gcf;
exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Helmholtz1D-%s-FreErr-Jacobi.pdf',index_num,mtd))
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Magnitude of ModeWise Error','Renderer', 'painters', 'Position', [0 0 700 350]);

iter_index=repmat((1:1:iters+1)',size(FreErr_Hybrid,2),1);
mode_index=reshape(repmat(1:1:size(FreErr_Hybrid,2),size(FreErr_Hybrid,1),1),size(FreErr_Hybrid,1)*size(FreErr_Hybrid,2),1);
Ft = TriScatteredInterp(iter_index,mode_index,reshape(FreErr_Hybrid,size(FreErr_Hybrid,1)*size(FreErr_Hybrid,2),1));
[qx,qy] = meshgrid(1:0.1:iters+1,1:0.1:size(FreErr_Hybrid,2));
qz = Ft(qx,qy);
x=[1 iters+1];
y=[1 size(FreErr_Hybrid,2)];
imagesc(x,y,qz)

colormap winter
cb = colorbar();
cb.Location = "eastoutside";
cb.AxisLocation = "out";
cb.Ruler.Exponent = 0;
set(gca,'ColorScale','log')
set(cb,'TickLabelInterpreter','latex')
axis xy
axis off
set(gca,'FontSize',18);
set(gcf,'color','w');
hold on
caxis([min_E max_E])

axis on
xticks(1:10:51)
xticklabels({'1','11','21','31','41','51'})
xlabel('index $k$ of iterations','interpreter','latex');
ylabel('index $j$ of frequencies','interpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ax = gcf;
exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Helmholtz1D-%s-FreErr-Hybrid.pdf',index_num,mtd));
%---------------------%

end

