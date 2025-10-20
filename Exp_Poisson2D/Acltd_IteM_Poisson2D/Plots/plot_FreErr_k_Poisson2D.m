function plot_FreErr_k_Poisson2D(FreErr_Hybrid,FreErr_Jacobi,k,filename,mtd)

    %---------------------%
    h=figure('NumberTitle','off','Name',sprintf('ModeErr during Iterations k=%d',k),'Renderer', 'painters', 'Position', [0 0 600 350]);
    semilogy(1:1:size(FreErr_Jacobi,2), FreErr_Jacobi(k,:), 'b-','LineWidth',1)
    hold on
    semilogy(1:1:size(FreErr_Hybrid,2), FreErr_Hybrid(k,:), 'r-','LineWidth',1 )
    axis xy
    xlim([1 size(FreErr_Jacobi,2)]);
    xticks(1:20:121)
    xticklabels({'1','21','41','61','81','101','121'})
    set(gca,'FontSize',18);
    set(gcf,'color','w');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    hl=legend(sprintf('${M}^{[k],%d}_{\\mathrm{Jacobi}}$',k),sprintf('${M}^{[k],%d}_{\\mathrm{Hybrid}}$',k),'interpreter','latex', 'orientation','horizontal','location', 'east');

    ax = gcf;
    exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-FreErr-k=%d.pdf',filename,mtd,k));
    %---------------------%

end