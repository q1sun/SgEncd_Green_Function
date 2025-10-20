function  plot_FreErr_k_Poisson1D(err_frq_Hybrid,err_frq_Jacobi,k,index_num,mtd)

    % plot the curves at same mode k
 	format short e
    number_of_iterations=size(err_frq_Hybrid,2)-1;

    %---------------------%
    h=figure('NumberTitle','off','Name',sprintf('ModeErr during Iterations k=%d',k),'Renderer', 'painters', 'Position', [0 0 600 350]);
    semilogy(1:1:number_of_iterations+1, err_frq_Jacobi(k,:), 'b-','LineWidth',1 )
    hold on
    semilogy(1:1:number_of_iterations+1, err_frq_Hybrid(k,:), 'r-','LineWidth',1)
    axis xy
    xlim([1 number_of_iterations+1]);
    xticks(1:10:51)
    xticklabels({'1','11','21','31','41','51'})
    set(gca,'FontSize',18);
    set(gcf,'color','w');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    hl=legend(sprintf('${M}^{[k],%d}_{\\mathrm{Jacobi}}$',k),sprintf('${M}^{[k],%d}_{\\mathrm{Hybrid}}$',k),'interpreter','latex', 'orientation','horizontal','location', 'east');
    hl.Position(2) = hl.Position(2) + 0.2;

    ax = gcf;
    exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Poisson1D-%s-FreErr-k=%d.pdf',index_num,mtd,k))
    %---------------------%
    
end

