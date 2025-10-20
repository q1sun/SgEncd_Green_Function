function plot_eigenpair_Green_Poisson1D( Mu_deep, Phi_deep, Mu_exact, Phi_exact, Mu_RelErr, V_basis, Phi_L2RelErr )

format short e


%% 1. eigenvalue
number_of_eigenvalues_begin = 1;
number_of_eigenvalues_end = 50; 

%---------------------%
%---------------------%
h=figure('NumberTitle','off','Name','Eigenvalues of Green Function','Renderer', 'painters', 'Position', [0 0 500 350]);
semilogy( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_exact(number_of_eigenvalues_begin:number_of_eigenvalues_end ), 'bo','LineWidth',0.5)%, 'MarkerSize',5)
hold on
semilogy( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_deep(number_of_eigenvalues_begin:number_of_eigenvalues_end ), 'r*','LineWidth',0.5)%, 'MarkerSize',5)

yticks([10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])
ylim([10^(-5) 10^(0)])
axis xy
set(gca,'FontSize',18);
set(gcf,'color','w')
set(gca,'XLim',[number_of_eigenvalues_begin number_of_eigenvalues_end]);
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('eigenvalue','interpreter','latex');
legend({'$\mu_j$', '$\check{\mu}_j$'}, 'location', 'northeast', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task3/fig-Poisson1D-eigenvalues-1to50.pdf')

%---------------------%
h=figure('NumberTitle','off','Name','Eigenvalues of Green Function (relative error)','Renderer', 'painters', 'Position', [0 0 500 350]);
plot( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_RelErr(number_of_eigenvalues_begin:number_of_eigenvalues_end), '-+','MarkerSize',4,'color', [0.4660 0.6740 0.1880]);
axis xy
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('relative error','interpreter','latex');
set(gca,'FontSize',18);
set(gca,'XLim',[number_of_eigenvalues_begin number_of_eigenvalues_end,]);
set(gcf,'color','w')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
legend({'$\frac{|\mu_j - \check{\mu}_j|}{|\mu_j|}$'}, 'location', 'southeast', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

ax = gcf;
exportgraphics(ax, 'Figures/Task3/fig-Poisson1D-eigenvalues-RelErr-1to50.pdf')
%---------------------%
%---------------------%

number_of_eigenvalues_begin = 150; 
number_of_eigenvalues_end = 200; 

%---------------------%
%---------------------%
h=figure('NumberTitle','off','Name','Eigenvalues of Green Function','Renderer', 'painters', 'Position', [0 0 500 350]);
semilogy( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_exact(number_of_eigenvalues_begin:number_of_eigenvalues_end ), 'bo','LineWidth',0.5)%, 'MarkerSize',5)
hold on
semilogy( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_deep(number_of_eigenvalues_begin:number_of_eigenvalues_end ), 'r*','LineWidth',0.5)%, 'MarkerSize',5)

axis xy
set(gca,'FontSize',18);
set(gcf,'color','w')
set(gca,'XLim',[number_of_eigenvalues_begin number_of_eigenvalues_end]);
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('eigenvalue','interpreter','latex');
legend({'$\mu_j$', '$\check{\mu}_j$'}, 'location', 'northeast', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = gcf;
exportgraphics(ax, 'Figures/Task3/fig-Poisson1D-eigenvalues-151to200.pdf')

%---------------------%
h=figure('NumberTitle','off','Name','Eigenvalues of Green Function (relative error)','Renderer', 'painters', 'Position', [0 0 500 350]);
plot( number_of_eigenvalues_begin:number_of_eigenvalues_end, Mu_RelErr(number_of_eigenvalues_begin:number_of_eigenvalues_end), '-+','MarkerSize',4,'color', [0.4660 0.6740 0.1880]);
axis xy
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('relative error','interpreter','latex');
set(gca,'FontSize',18);
set(gca,'XLim',[number_of_eigenvalues_begin number_of_eigenvalues_end,]);
set(gcf,'color','w')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
legend({'$\frac{|\mu_j - \check{\mu}_j|}{|\mu_j|}$'}, 'location', 'southeast', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

ax = gcf;
exportgraphics(ax, 'Figures/Task3/fig-Poisson1D-eigenvalues-RelErr-151to200.pdf')
%---------------------%
%---------------------%



%% 2. eigenfunction
Phi_exact = Phi_exact * sqrt(2);
Phi_deep = Phi_deep * sqrt(2);

%---------------------%
h=figure('NumberTitle','off','Name','First 5 eigenfunctions of deep Green Function','Renderer', 'painters', 'Position', [0 0 500 350]);
p1 = plot( V_basis(:), sign(Phi_deep(10,1)) * Phi_deep(:,1), '--','color',[1 0 0 1], 'LineWidth', 2 );
hold on
p5 = plot( V_basis(:), sign(Phi_exact(10,1)) * Phi_exact(:,1), 'color', [1, 0.41176, 0.70588, 0.5], 'LineWidth', 2 );
hold on
p2 = plot( V_basis(:), sign(Phi_deep(10,3)) * Phi_deep(:,3), '--','color',[0 0.39216 0 1], 'LineWidth', 2 );
hold on
p6 = plot( V_basis(:), sign(Phi_exact(10,3)) * Phi_exact(:,3), 'color', [0, 1, 0, 0.5], 'LineWidth', 2 );
hold on
p3 = plot( V_basis(:), sign(Phi_deep(10,9)) * Phi_deep(:,9), '--', 'color', [0 0 1 1], 'LineWidth', 2 );
hold on
p7 = plot( V_basis(:), sign(Phi_exact(10,9)) * Phi_exact(:,9), 'color', [0.3010 0.7450 0.9330 0.5], 'LineWidth', 2 );
p1.Color(4) = 1;  p2.Color(4) = 1;  p3.Color(4) = 1; 
axis xy
xlabel('$x$','interpreter','latex');
ylabel('eigenfunction','interpreter','latex');
set(gca,'FontSize',18);
set(gcf,'color','w')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ylim([-sqrt(2) sqrt(2)])
legend({'$\check{\phi}_1(x)$', '$\phi_1(x)$', '$\check{\phi}_3(x)$', '$\phi_3(x)$', '$\check{\phi}_{9}(x)$',  ...
    '$\phi_{9}(x)$'}, 'NumColumns', 3, 'location', 'south', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

ax = gcf;
exportgraphics(ax, 'Figures/Task3/fig-Poisson1D-eigenfunctions-deep-exact.pdf')
%---------------------%

% relative error for eigenfunction

number_of_eigenvalues_begin = 1; % show eigenvalues
number_of_eigenvalues_end = 200; % show eigenvalues

%---------------------%
h=figure('NumberTitle','off','Name','L2 Error of Eigenfunctions (relative error)','Renderer', 'painters', 'Position', [0 0 570 360]);
semilogy( number_of_eigenvalues_begin:number_of_eigenvalues_end, Phi_L2RelErr(number_of_eigenvalues_begin:number_of_eigenvalues_end), '-+','MarkerSize',4,'color',[0.44 0.5 0.41]);
axis xy
xlabel('index $k$ of eigenfunction','interpreter','latex');
ylabel('relative error','interpreter','latex');
set(gca,'FontSize',18);
set(gca,'XLim',[number_of_eigenvalues_begin number_of_eigenvalues_end]);
set(gca,'FontSize',18);
set(gcf,'color','w')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
legend({'$\frac{\| \phi_j - \check{\phi}_j \|_{L^2(\Omega)}}{\|\phi_j\|_{L^2(\Omega)}}$'}, 'location', 'southeast', 'interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

ax = gcf;
saveas(ax, 'Figures/Task3/fig-Poisson1D-eigenfunctions-RelErr-1to200.pdf')
%---------------------%

end

