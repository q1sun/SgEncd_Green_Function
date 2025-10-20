function run_DownStreamTask_SgEncd_Green_Poisson1D

format short e

% task-1: plot SgEncd Green's function
% task-2: compute solution using SgEncd Green's function
% task-3: recover eigenpairs using SgEncd Green's function

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
left = 0; right = 1; % domain
Green_exact = @(x,y) x .*(1-y) .* (x<=y) + y .* (1-x) .* (x>y); % Green's function

u_exact = @(x) 10 * x - 10 * x.^2 + 0.5 * sin( 20 * pi * x.^3 ); % test example 
f_exact = @(x) - 60 * pi * x .* cos(20 * pi * x.^3) + 1800 * pi^2 * x.^4 .* sin( 20 * pi * x.^3) + 20; % test example 

phi_exact = @(x,k) sin( k .* pi .* x);
%-------------------------------------------------------%
% 0-2. subroutines
%-------------------------------------------------------%
addpath('Data/');
addpath('Plots/');
addpath('Task1_Plot_SgEncdGreen/');
addpath('Task2_Represent_Solution/');
addpath('Task3_Eigenvalue_Problem/');
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
%-------------------------------------------------------%
% task-3. recover eigenpairs using SgEncd Green's function
%-------------------------------------------------------%
h = 1/2^10; % mesh size
basis_type = 102; % type of finite elments
[P_mesh, T_mesh, V_basis, T_basis] = Generate_Information_Matrices_Linear_Quadratic_1D(left,right,h,basis_type);
task3_mesh_Gpts = assemble_input_meshgrid(V_basis, V_basis);

save('Data/task3_mesh_Gpts.mat', 'task3_mesh_Gpts'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run python to load SgEncd Green's function values

%% 2. load SgEncd Green's function for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% task-1 . plot SgEncd Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_FigMesh_Value.mat'); 
% plot exact and predicted Green's function on xy full domain
exact_G_meshXY = Green_exact( task1_mesh_Gpts(1,:), task1_mesh_Gpts(2,:) );
plot_Green_function_Poisson1D(exact_G_meshXY, NN_G_meshXY', task1_mesh_Gpts);
% plot exact and predicted Green's function at y = 0.3
exact_G_fixedY_trace = Green_exact( mesh_x, 0.3 );
plot_fixedY_Green_function_Poisson1D(exact_G_fixedY_trace, NN_G_fixedY', NN_G_fixedY_trace', task1_mesh_Gpts, mesh_x);
%-------------------------------------------------------%
% task-2. compute solution using deep Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_Qpts_Value.mat')
f_Qpts = f_exact( task2_mesh_Qpts(2,:));
temp = NN_G_Qpts .* f_Qpts' .* repmat(W_Qpts, number_of_elements + 1, 1);
NN_u_meshpts = sum( reshape(temp, number_of_elements * number_of_quadrature_points, number_of_elements+1) );
exact_u_meshpts = u_exact( task2_V_mesh(:) );
f_meshpts_exact = f_exact( task2_V_mesh(:) );

plot_solution_representation_Poisson1D( task2_V_mesh, exact_u_meshpts, NN_u_meshpts, f_meshpts_exact );
%-------------------------------------------------------%
% task-3. recover eigenpairs using deep Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_Fem_Value.mat')
% compute eigenpairs of the SgEncd Green's function
number_of_eigs = 200;
NN_G_Fem = reshape( NN_G_Fem, size(V_basis,2), size(V_basis,2));
[Mu_neural, Phi_neural] = eigensolver_Green_Poisson1D(NN_G_Fem, basis_type, P_mesh, T_mesh, V_basis, T_basis, number_of_eigs);
% obtain exact eigenpairs
Mu_exact = ( (1:number_of_eigs).^2 * pi^2 ).^(-1);
for i = 1 : number_of_eigs
    Phi_exact(:, i) = phi_exact( V_basis(:), i );
end
% relative error of eigenpairs
Mu_RelErr = abs(Mu_exact' - Mu_neural(1:number_of_eigs)) ./ Mu_exact';
Phi_L2RelErr = compute_eigenfunction_L2Err_Poisson1D(Phi_neural, number_of_eigs, P_mesh, T_mesh, T_basis, basis_type);

% plot eigenpairs
plot_eigenpair_Green_Poisson1D( Mu_neural, Phi_neural, Mu_exact, Phi_exact, Mu_RelErr, V_basis, Phi_L2RelErr );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







end

