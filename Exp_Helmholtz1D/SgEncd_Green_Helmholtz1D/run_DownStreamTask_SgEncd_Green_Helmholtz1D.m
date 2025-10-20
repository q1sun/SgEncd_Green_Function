function run_DownStreamTask_SgEncd_Green_Helmholtz1D

format short e

% task-1: plot SgEncd Green's function
% task-2: compute solution using SgEncd Green's function

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
left = 0; right = 1; % domain
% base functions
k_exact = @(x) 15 * sin( 10 * x );
u_exact = @(x) 10 * x - 10 * x.^2 + 0.5 * sin( 20 * pi * x.^3 ); 
% first and second derivatives: du/dx, d^2u/d^2x
du_exact = @(x) 10 - 20 * x + 30 * pi * x.^2 .* cos(20 * pi * x.^3);
d2u_exact = @(x) -20 + 60 * pi * x .* cos(20 * pi * x.^ 3) - 1800 * pi^2 * x.^4 .* sin(20 * pi * x .^3);
% forcing term f for the equation: -((x-2)^2 * u'') - (2*(x-2) * u') - (k^2 * u)
f_exact = @(x) -2 * (x-2) .* du_exact(x) - (x-2) .^ 2 .* d2u_exact(x) - k_exact(x) .^2 .* u_exact(x); 
%-------------------------------------------------------%
% 0-2. subroutines
%-------------------------------------------------------%
addpath('Data/');
addpath('Plots/');
addpath('Task1_Plot_SgEncdGreen/');
addpath('Task2_Represent_Solution/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Generate InputData of SgEncd Green's function for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% task-1 . plot SgEncd Green's function
%-------------------------------------------------------%
num_pts = 200; % number of mesh points for figure
mesh_x = linspace(left, right, num_pts);
mesh_y = linspace(left, right, num_pts);
task1_mesh_Gpts = assemble_input_meshgrid(mesh_x, mesh_y);

save('Data/task1_mesh_Gpts.mat', 'mesh_x', 'mesh_y', 'task1_mesh_Gpts');
%-------------------------------------------------------%
% task-2. compute solution using SgEncd Green's function
%-------------------------------------------------------%
number_of_elements = 2^10; % piecewise quadrature rule
number_of_quadrature_points = 8;
[task2_V_mesh, V_Qpts, W_Qpts] = generate_mesh_quadrature_points(left, right, number_of_elements, number_of_quadrature_points);
task2_mesh_Qpts = assemble_input_meshgrid(task2_V_mesh, V_Qpts);

save('Data/task2_mesh_Qpts.mat', 'task2_mesh_Qpts');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run python to load SgEncd Green's function values

%% 2. load SgEncd Green's function for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% task-1 . plot SgEncd Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_FigMesh_Value.mat'); 
% plot predicted Green's function on xy full domain / for fixed y = 0.3
plot_Green_function_Helmholtz1D(NN_G_meshXY', task1_mesh_Gpts);
plot_fixedY_Green_function_Helmholtz1D(NN_G_fixedY', NN_G_fixedY_trace', task1_mesh_Gpts, mesh_x);
%-------------------------------------------------------%
% task-2. compute solution using deep Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_Qpts_Value.mat')
f_Qpts = f_exact( task2_mesh_Qpts(2,:));
temp = NN_G_Qpts .* f_Qpts' .* repmat(W_Qpts, number_of_elements + 1, 1);
NN_u_meshpts = sum( reshape(temp, number_of_elements * number_of_quadrature_points, number_of_elements+1) );
exact_u_meshpts = u_exact( task2_V_mesh(:) );
f_meshpts_exact = f_exact( task2_V_mesh(:) );

plot_solution_representation_Helmholtz1D( task2_V_mesh, exact_u_meshpts, NN_u_meshpts, f_meshpts_exact );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







end



