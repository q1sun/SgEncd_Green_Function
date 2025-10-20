function plot_eigenpair_Green_Poisson2D(lambda_pred, lambda_exact, num_eigvals)

%% plot eigenvalues
num_eigvals_first = 50; % number of eigenvalues being plotted
%---------------------%
figure('NumberTitle','off','Name','eigenvalues-deep', 'Position', [0 0 500 350]);

semilogy(1:1:num_eigvals_first,lambda_exact(1:num_eigvals_first),'bo','LineWidth',0.5);
hold on 
semilogy(1:1:num_eigvals_first,lambda_pred(1:num_eigvals_first),'r*','LineWidth',0.5);

yticks([10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])
ylim([10^(-5) 10^(0)])
axis xy
set(gca,'FontSize',18);
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('eigenvalue','Interpreter','latex');
legend({'$\mu_j$', '$\check{\mu}_j$'}, 'location', 'northeast', 'interpreter', 'latex');
xlim([1 num_eigvals_first])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(gcf,'color','w')
ax = gcf;   
exportgraphics(ax, 'Figures/Task3/fig-Poisson2D-eigenvalues-1to50.pdf')
%---------------------%

%---------------------%
figure('NumberTitle','off','Name','eigenvalues-deep', 'Position', [0 0 500 350]);

semilogy(251:1:300,lambda_exact(251:300),'bo','LineWidth',0.5);
hold on 
semilogy(251:1:300,lambda_pred(251:300),'r*','LineWidth',0.5);

yticks([10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])
ylim([10^(-5) 10^(0)])
axis xy
set(gca,'FontSize',18);
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('eigenvalue','Interpreter','latex');
legend({'$\mu_j$', '$\check{\mu}_j$'}, 'location', 'northeast', 'interpreter', 'latex');
xlim([251 300])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(gcf,'color','w')
ax = gcf;   
exportgraphics(ax, 'Figures/Task3/fig-Poisson2D-eigenvalues-251to300.pdf')
%---------------------%

%---------------------%
figure('NumberTitle','off','Name','eigenvalues-error', 'Position', [0 0 500 350]);

lambda_err = lambda_pred(1:num_eigvals)-lambda_exact(1:num_eigvals);
lambda_relative_err = abs(lambda_err) ./ lambda_exact(1:num_eigvals);
plot(1:num_eigvals, lambda_relative_err,'-+','MarkerSize',4,'color', [0.4660 0.6740 0.1880]);
axis xy
xlabel('index $j$ of eigenvalue','interpreter','latex');
ylabel('relative error','Interpreter','latex');
legend({'$\frac{|\mu_j - \check{\mu}_j|}{|\mu_j|}$'}, 'location', 'north', 'interpreter', 'latex');
xlim([1 num_eigvals])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(gcf,'color','w')
set(gca,'FontSize',18);
ax = gcf;

exportgraphics(ax, 'Figures/Task3/fig-Poisson2D-eigenvalues-RelErr-1to300.pdf')
%---------------------%

end