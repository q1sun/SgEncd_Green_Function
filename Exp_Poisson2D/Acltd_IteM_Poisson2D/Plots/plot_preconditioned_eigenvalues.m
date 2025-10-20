function plot_preconditioned_eigenvalues(Ah, NN_Preconditioner,index_grid)

iterator_Green = NN_Preconditioner * Ah;
[eigVec_Green, eigVal_Prectd] = eig(iterator_Green);
eigVal_Prectd = diag(eigVal_Prectd);
num_condition_prectd = max(eigVal_Prectd) / min(eigVal_Prectd);

[eigVec_A, eigVal_A] = eig(Ah);
eigVal_A = diag(eigVal_A);
num_condition_A = max(eigVal_A) / min(eigVal_A);

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot(eigVal_Prectd',0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0,5]);
ylim([0,1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,1,2,3,4,5]);
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, sprintf('Figure/Meshsize_h/fig-Poisson2D-eigenval-green-h%d.pdf',index_grid));
%---------------------%

fprintf('h%d: condition number(A) = %f \n',index_grid, num_condition_A);
fprintf('h%d: condition number = %f \n',index_grid, num_condition_prectd);

