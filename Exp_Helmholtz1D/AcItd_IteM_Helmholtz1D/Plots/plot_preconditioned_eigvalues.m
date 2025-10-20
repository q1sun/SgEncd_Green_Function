function plot_preconditioned_eigvalues(prectd_matrix_h6, prectd_matrix_h8, prectd_matrix_h10, prectd_matrix_h12)

% input: the preconditioned matrix BA of different sizes

eigVal_prectd_h6 = eig(prectd_matrix_h6);
kappa_h6 = max(eigVal_prectd_h6) / min(eigVal_prectd_h6);

eigVal_prectd_h8= eig(prectd_matrix_h8);
kappa_h8 = max(eigVal_prectd_h8) / min(eigVal_prectd_h8);

eigVal_prectd_h10 = eig(prectd_matrix_h10);
kappa_h10 = max(eigVal_prectd_h10) / min(eigVal_prectd_h10);

eigVal_prectd_h12 = eig(prectd_matrix_h12);
kappa_h12 = max(eigVal_prectd_h12) / min(eigVal_prectd_h12);

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot(eigVal_prectd_h6',0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0,5]);
ylim([0,0.1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,1,3,5]);
set(gca,'XTickLabel',{'$0$','$1$','$3$', '$5$'});
% set(gca,'xtick',[1e2,1e3,1e4]);
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, 'Figure/Meshsize_h/fig-Helmholtz1D-eigenval-green-h6.pdf');
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot(eigVal_prectd_h8',0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0,5]);
ylim([0,0.1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,1,3,5]);
set(gca,'XTickLabel',{'$0$','$1$','$3$', '$5$'});
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, 'Figure/Meshsize_h/fig-Helmholtz1D-eigenval-green-h8.pdf');
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot(eigVal_prectd_h10',0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0,5]);
ylim([0,0.1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,1,3,5]);
set(gca,'XTickLabel',{'$0$','$1$','$3$', '$5$'});
% set(gca,'xtick',[1e2,1e3,1e4]);
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, 'Figure/Meshsize_h/fig-Helmholtz1D-eigenval-green-h10.pdf');
%---------------------%

%---------------------%
h=figure('NumberTitle','off','Name','Solution','Renderer', 'painters', 'Position', [0 0 500 50]);
%crosses
plot(eigVal_prectd_h12',0,'ro','markersize',4)
yline(0,'linewidth',1.5)
xlim([0,5]);
ylim([0,0.1]);
set(gca,'ytick',[]);
% delete the axis
box off;
ax=gca;
pause(1e-6)
ax.XRuler.Axle.Visible='off';
ax.YRuler.Axle.Visible='off';
ax.TickDir = 'none';
set(gca,'xtick',[0,1,3,5]);
set(gca,'XTickLabel',{'$0$','$1$','$3$', '$5$'});
set(gca,'FontSize',18);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
ax = gcf;
exportgraphics(ax, 'Figure/Meshsize_h/fig-Helmholtz1D-eigenval-green-h12.pdf');
%---------------------%

fprintf('h=1/2^6: condition number = %f \n',kappa_h6);
fprintf('h=1/2^8: condition number = %f \n',kappa_h8);
fprintf('h=1/2^10: condition number = %f \n',kappa_h10);
fprintf('h=1/2^12: condition number = %f \n',kappa_h12);





