function  plot_IterErr_MG_Helmholtz1D(eh_L2_2gMG, rh_L2_2gMG,...
    eh_L2_5gMG, rh_L2_5gMG, eh_L2_HbMG, rh_L2_HbMG,...
    index_num,mtd)
    
    % plot the error curves of the two different methods
    format short e
    max_iter=size(eh_L2_5gMG,2)-1;
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Error for Solution Prediction','Renderer', 'painters', 'Position', [0 0 700 350]);
    semilogy(1:1:max_iter+1, eh_L2_HbMG, 'LineWidth',1.5 );
    hold on
    semilogy(1:1:max_iter+1, eh_L2_2gMG,  'LineWidth',1.5 );
    hold on
    semilogy(1:1:max_iter+1, eh_L2_5gMG, 'm', 'LineWidth',1.5 );
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    lgd = legend('2-grid MG',...
        '5-grid MG',...
        '2-grid Hybrid MG',....
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    lgd.Position(2) = lgd.Position(2) + 0.1;
    xlim([1 21])
    xticks(1:5:21)
    xticklabels({'1','6','11','16','21'})
    xlabel('number of cycles','interpreter','latex');
    ylabel('solution error $\| E^{[k]} \|_2$','interpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Helmholtz1D-%s-IteErr-2Norm.pdf',index_num,mtd));
    %---------------------%
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Residual for Solution Prediction','Renderer', 'painters', 'Position', [0 0 700 350]);
    semilogy(1:1:max_iter+1, rh_L2_HbMG, 'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_2gMG,  'LineWidth',1.5 )
    hold on
    semilogy(1:1:max_iter+1, rh_L2_5gMG,'m', 'LineWidth',1.5 )
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w')
    legend('2-grid MG',...
        '5-grid MG',...
        '2-grid Hybrid MG',....
        'interpreter', 'latex', 'location', 'east','Fontsize',15);
    xlim([1 21])
    xticks(1:5:21)
    xticklabels({'1','6','11','16','21'})
    xlabel('number of cycles','interpreter','latex');
    ylabel('residual error $\| R^{[k]} \|_2$','interpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Helmholtz1D-%s-ResErr-2Norm.pdf',index_num,mtd));
    %---------------------%
    
end


