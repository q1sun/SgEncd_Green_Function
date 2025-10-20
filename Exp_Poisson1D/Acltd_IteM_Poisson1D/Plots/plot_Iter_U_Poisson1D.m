function plot_Iter_U_Poisson_line1D(V_x_Gpts,uh_Hybrid,uh_Jacobi,uh,eh_Hybrid,eh_Jacobi,index_num,mtd)

    format long e

    % Attention: Here Xh,Vh,V_exact includes the boundary points.
    number_of_iteration=size(uh_Hybrid,2);
    uh_Hybrid=[zeros(1,number_of_iteration);uh_Hybrid;zeros(1,number_of_iteration)];
    uh_Jacobi=[zeros(1,number_of_iteration);uh_Jacobi;zeros(1,number_of_iteration)];
    eh_Hybrid=[zeros(1,number_of_iteration);eh_Hybrid;zeros(1,number_of_iteration)];
    eh_Jacobi=[zeros(1,number_of_iteration);eh_Jacobi;zeros(1,number_of_iteration)];
    uh=[0;uh;0];
    
    %---------------------%
    h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 420]);
    plot(V_x_Gpts, uh, 'k-','LineWidth',1.5 )
    hold on
    plot(V_x_Gpts,uh_Jacobi(:,1),'--','color',[0.5 0.5 0.5],'LineWidth',1.5)
    hold on
    plot(V_x_Gpts, uh_Hybrid(:,end), 'r--','LineWidth',1.2 )
    hold on
    plot(V_x_Gpts, uh_Jacobi(:,end), 'b--','LineWidth',1.1 )
    axis xy
    set(gca,'FontSize',18);
    set(gcf,'color','w');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    xl = xlabel('$x$','interpreter','latex');
    xl.Position(2) = xl.Position(2) +  0.3;
    ax = gca;
    ax.XAxis.TickLabelGapMultiplier = 0.1; 
    ylabel('solution value at $k$-th iteration','interpreter','latex');
    h1=legend({'$U$','$U^{[1]}$','$U_{\mathrm{Jacobi}}^{[51]}$','$U_{\mathrm{Hybrid}}^{[51]}$'}, 'location', 'south', 'interpreter', 'latex','FontSize',14);
    set(h1,'Orientation','horizon');
    set(h1,'position',[0.159434815062664,0.171156893819335,0.714419317943913,0.081622926255602]);
    h1.Position(2) = h1.Position(2) -  0.05;
    xticks([0 1])
    xticklabels({'0','1'})
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Poisson1D-%s-U-k25.pdf',index_num,mtd));
    %---------------------%

    %---------------------%
    h=figure('NumberTitle','off','Name','Hybrid-Process of Nodewise Error','Renderer', 'painters', 'Position', [0 0 750 350]);
    t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
    nexttile
    plot(V_x_Gpts, eh_Jacobi(:,1), 'b-','LineWidth',0.8);
    l1=legend('$E_{\mathrm{Jacobi}}^{[1]}$','FontSize',12,'location','South','interpreter', 'latex');
    set(gca,'FontSize',14);
    xticks([0 1])
    xticklabels({'0','1'})
    nexttile
    plot(V_x_Gpts, eh_Hybrid(:,1), 'r-','LineWidth',0.8);
    l4 = legend('$E_{\mathrm{Hybrid}}^{[1]}$','FontSize',12,'location', 'South', 'interpreter', 'latex');
    set(gca,'FontSize',14);
    xticks([0 1])
    xticklabels({'0','1'})
    nexttile
    plot(V_x_Gpts, eh_Jacobi(:,20), 'b-','LineWidth',0.8);
    l2=legend('$E_{\mathrm{Jacobi}}^{[21]}$','FontSize',12,'location', 'South', 'interpreter', 'latex');
    set(gca,'FontSize',14);
    xticks([0 1])
    xticklabels({'0','1'})
    nexttile
    plot(V_x_Gpts, eh_Hybrid(:,20), 'r-','LineWidth',0.8);
    l5 = legend('$E_{\mathrm{Hybrid}}^{[21]}$','FontSize',12,'location', 'South', 'interpreter', 'latex');
    set(gca,'FontSize',14);
    xticks([0 1])
    xticklabels({'0','1'})
    nexttile
    plot(V_x_Gpts, eh_Jacobi(:,30), 'b-','LineWidth',0.8 );
    set(gca,'FontSize',14);
    xticks([0 1])
    xticklabels({'0','1'})
    l3 = legend('$E_{\mathrm{Jacobi}}^{[31]}$','FontSize',12,'location', 'south', 'interpreter', 'latex');
    l1.Position(1) = l1.Position(1) +  0.02;
    l1.Position(2) = l1.Position(2) -  0.02;
    l2.Position(2) = l2.Position(2) -  0.02;
    l2.Position(1) = l2.Position(1) +  0.02;
    l3.Position(2) = l3.Position(2) -  0.02;
    l3.Position(1) = l3.Position(1) +  0.04;
    nexttile
    plot(V_x_Gpts, eh_Hybrid(:,30), 'r-','LineWidth',0.8 );
    l6 = legend('$E_{\mathrm{Hybrid}}^{[31]}$','FontSize',12,'location', 'South', 'interpreter', 'latex');
    set(gca,'FontSize',14);
    l4.Position(2) = l4.Position(2) -  0.02;
    l5.Position(2) = l5.Position(2) -  0.02;
    l6.Position(2) = l6.Position(2) -  0.02;
    axis xy
    xlabel(t,'','interpreter','latex','FontSize',18);
    ylabel(t,'solution error at $k$-th iteration','interpreter','latex','FontSize',18)
    annotation('textbox', [0.24, 0.00, 0.1, 0.03], ...
           'String', '$x$', ...
           'Interpreter', 'latex', ...
           'FontSize', 18, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', ...
           'EdgeColor', 'none'); 
    annotation('textbox', [0.72, 0.00, 0.1, 0.03], ...
           'String', '$x$', ...
           'Interpreter', 'latex', ...
           'FontSize', 18, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', ...
           'EdgeColor', 'none'); 
    set(gcf,'color','w')
    xticks([0 1])
    xticklabels({'0','1'})
    set(gca,'FontSize',14);
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    ax = gcf;
    exportgraphics(ax, sprintf('Figure/Meshsize_2e-%d/fig-Poisson1D-%s-U-Smooth-IteErr.pdf',index_num,mtd));
    %---------------------%
    
    
end



