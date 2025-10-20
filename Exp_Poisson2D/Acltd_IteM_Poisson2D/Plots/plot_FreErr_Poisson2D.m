function plot_FreErr_Poisson2D(FreErr_Hybrid,FreErr_Jacobi,filename,mtd)

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
x=[1 size(FreErr_Jacobi,1)+1];
y=[1 size(FreErr_Jacobi,2)];
%qz=Modewise_Err_Jacobi';
imagesc(x,y,qz)
%axis([1 size(EH_Jacobi,1) 1 size(EH_Jacobi,2)]);

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
xticks(1:20:121)
xticklabels({'1','21','41','61','81','101','121'})
xlabel('index $k$ of iterations','interpreter','latex');
ylabel('index $j$ of frequencies','interpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ax = gcf;
exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-FreErr-Jacobi.pdf',filename,mtd))
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
%qz=Modewise_Err_Hybrid';
imagesc(x,y,qz)
%axis([1 size(Modewise_Err_Hybrid,1) 1 size(Modewise_Err_Hybrid,2)]);

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
xticks(1:20:121)
xticklabels({'1','21','41','61','81','101','121'})
xlabel('index $k$ of iterations','interpreter','latex');
ylabel('index $j$ of frequencies','interpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ax = gcf;
exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-FreErr-Hybrid.pdf',filename,mtd));
%---------------------%

end