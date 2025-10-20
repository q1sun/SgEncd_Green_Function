function  plot_IterErr_MG_Poisson2D(eh_L2_Jacobi, rh_L2_Jacobi,...
    eh_L2_Hybrid_2L, rh_L2_Hybrid_2L, eh_L2_Jacobi_2L, rh_L2_Jacobi_2L,...
    filename,mtd)
    
    % plot the error curves of the two different methods
    format short e
    max_iter=size(eh_L2_Hybrid_2L,2)-1;
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Error for Solution Prediction','Renderer', 'painters', 'Position', [0 0 700 350]);
    p1 = semilogy(1:1:max_iter+1, eh_L2_Jacobi_2L, 'LineWidth',1.5 );
    hold on
    p1 = semilogy(1:1:max_iter+1, eh_L2_Jacobi,  'LineWidth',1.5 );
    hold on
    
    p2 = semilogy(1:1:max_iter+1, eh_L2_Hybrid_2L, 'm', 'LineWidth',1.5 );
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    lgd = legend('2-grid MG',...
        '4-grid MG',...
        '2-grid Hybrid MG',....
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    lgd.Position(2) = lgd.Position(2) + 0.1;
    xlim([1 31])
    xticks(1:5:31)
    xticklabels({'1','6','11','16','21','26','31'})
    xlabel('number of cycles','interpreter','latex');
    ylabel('solution error $\| E^{[k]} \|_2$','interpreter','latex');
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
    semilogy(1:1:max_iter+1, rh_L2_Jacobi_2L, 'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_Jacobi,  'LineWidth',1.5 )
    hold on
    
    semilogy(1:1:max_iter+1, rh_L2_Hybrid_2L,'m', 'LineWidth',1.5 )
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    lgd = legend('2-grid MG',...
        '4-grid MG',...
        '2-grid Hybrid MG',....
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    % lgd.Position(2) = lgd.Position(2) - 0.1;
    xlim([1 31])
    xticks(1:5:31)
    xticklabels({'1','6','11','16','21','26','31'})
    xlabel('number of cycles','interpreter','latex');
    ylabel('residual error $\| R^{[k]} \|_2$','interpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/%s/fig-Poisson2D-%s-ResErr-2Norm.pdf',filename,mtd));
    %---------------------%
    
end
