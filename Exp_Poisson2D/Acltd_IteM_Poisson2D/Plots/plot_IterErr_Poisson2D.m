function  plot_IterErr_Poisson2D(eh_L2_Hybrid, rh_L2_Hybrid,eh_L2_Jacobi, rh_L2_Jacobi,...
    eh_L2_Hybrid_2, rh_L2_Hybrid_2,...
    eh_L2_Hybrid_3, rh_L2_Hybrid_3,...
    eh_L2_Hybrid_4, rh_L2_Hybrid_4,...
    filename,mtd)
    
    % plot the error curves of the two different methods
    format short e
    max_iter=size(eh_L2_Hybrid,2)-1;
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Error for Solution Prediction','Renderer', 'painters', 'Position', [0 0 700 350]);
    p2 = semilogy(1:1:max_iter+1, eh_L2_Hybrid, 'r-','LineWidth',1.5 );
    hold on
    p3 = semilogy(1:1:max_iter+1, eh_L2_Hybrid_2,'g-', 'LineWidth',1.5 );
    hold on
    p4 = semilogy(1:1:max_iter+1, eh_L2_Hybrid_3, 'LineWidth',1.5 );
    hold on
    p5 = semilogy(1:1:max_iter+1, eh_L2_Hybrid_4, 'LineWidth',1.5 );
    hold on
    p1 = semilogy(1:1:max_iter+1, eh_L2_Jacobi, 'b-','LineWidth',1.5 );
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    lgd = legend('$K=2$',...
        '$K=5$',...
        '$K=10$',...
        '$K=20$',...
        '$\|E_{\mathrm{Jacobi}}^{[k]}\|_2$',...
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    lgd.Position(2) = lgd.Position(2) + 0.1;
    xlim([1 201])
    xticks(1:50:201)
    xticklabels({'1','51','101','151','201'})
    xlabel('index $k$ of iterations','interpreter','latex');
    ylabel('solution error $\| E_{\mathrm{Hybrid}}^{[k]} \|_2$','interpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-IteErr-2Norm.pdf',filename,mtd));
    %---------------------%
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Residual for Solution Prediction','Renderer', 'painters', 'Position', [0 0 700 350]);
    semilogy(1:1:max_iter+1, rh_L2_Hybrid, 'r-','LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_Hybrid_2,'g-', 'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_Hybrid_3, 'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_Hybrid_4, 'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_Jacobi, 'b-','LineWidth',1.5 )
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    lgd = legend('$K=2$',...
        '$K=5$',...
        '$K=10$',...
        '$K=20$',...
        '$\|R_{\mathrm{Jacobi}}^{[k]}\|_2$',...
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    lgd.Position(2) = lgd.Position(2) + 0.1;
    xlim([1 201])
    xticks(1:50:201)
    xticklabels({'1','51','101','151','201'})
    xlabel('index $k$ of iterations','interpreter','latex');
    ylabel('residual error $\| R_{\mathrm{Hybrid}}^{[k]} \|_2$','interpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-ResErr-2Norm.pdf',filename,mtd));
    %---------------------%
    
end
