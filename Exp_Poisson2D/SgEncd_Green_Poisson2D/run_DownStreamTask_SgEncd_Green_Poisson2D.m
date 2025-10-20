function run_DownStreamTask_SgEncd_Green_Poisson2D

format short e

% task-1: plot SgEncd Green's function
% task-2: compute solution using SgEncd Green's function
% task-3: recover eigenpairs using SgEncd Green's function

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
Green_exact = @(x,y)-0.5/pi*log(sqrt(sum((x-y).^2,2))./(sqrt(sum((x-y./sum(y.^2,2)).^2,2)).*sqrt(sum(y.^2,2)))); % Green's function

u_exact = @(x,y) - exp( x.^2 + y.^2 - 1 ) + 1; % test example
f_exact = @(x,y) 4 * exp( x.^2 + y.^2 - 1 ) .* (x.^2 + y.^2 + 1); % test example 
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
load('Data/MeshData_Disc2D.mat');
task1_Gpts_FixedY1 = [node; repmat([0;0], 1, size(node,2))];
task1_Gpts_FixedY2 = [node; repmat([0.7;0.2], 1, size(node,2))];
save('Data/task1_mesh_Gpts.mat', 'task1_Gpts_FixedY1', 'task1_Gpts_FixedY2');
%-------------------------------------------------------%
% task-2. compute solution using SgEncd Green's function
%-------------------------------------------------------%
number_of_elements = size(elem,2);
number_of_quadrature_points = 7;
[task2_V_mesh, V_Qpts, W_Qpts] = generate_mesh_quadrature_points(node,elem+1,number_of_quadrature_points);
task2_mesh_Qpts = assemble_input_meshgrid(task2_V_mesh, V_Qpts');

save('Data/task2_mesh_Qpts.mat', 'task2_mesh_Qpts');
%-------------------------------------------------------%
% task-3. recover eigenpairs using SgEncd Green's function
%-------------------------------------------------------%
num_pts = size(node,2);
% generate sample points within the ball for numerical integration
N_y = 5; % 5 gauss points 1d
eps = 1e-2; % radius of each ball
[task3_diag_Qpts, weights_diag] = generate_diagonal_gauss_points(node', eps, N_y);
task3_mesh_Gpts = assemble_input_meshgrid(node, node);
save('Data/task3_mesh_Gpts.mat', 'task3_diag_Qpts', 'task3_mesh_Gpts');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run python to load SgEncd Green's function values

%% 2. load SgEncd Green's function for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% task-1 . plot SgEncd Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_FigMesh_Value.mat');
% plot exact and predicted Green's function at fixed y
exact_G_fixedY1 = Green_exact(transpose(task1_Gpts_FixedY1(3:4,:)), transpose(task1_Gpts_FixedY1(1:2,:)));
exact_G_fixedY2 = Green_exact(transpose(task1_Gpts_FixedY2(1:2,:)), transpose(task1_Gpts_FixedY2(3:4,:)));
plot_fixedY_Green_function_Poisson2D(exact_G_fixedY1, double(NN_G_fixedY1), exact_G_fixedY2, double(NN_G_fixedY2), node);
%-------------------------------------------------------%
% task-2. compute solution using SgEncd Green's function
%-------------------------------------------------------%
f_Qpts = f_exact(V_Qpts(:,1),V_Qpts(:,2));

temp = [];
for i=0:1:6
    load(sprintf('Data/SgEncdGreen_Qpts_Value_%d.mat',i));
    temp = [temp; NN_G_Qpts];
end
NN_G_Qpts = double(temp);

NN_u_meshpts = reshape(NN_G_Qpts,size(task2_V_mesh,2),number_of_elements * number_of_quadrature_points)* (f_Qpts.* W_Qpts);
exact_u_meshpts = u_exact( task2_V_mesh(1,:)',task2_V_mesh(2,:)');
plot_solution_representation_Poisson2D( task2_V_mesh, exact_u_meshpts, NN_u_meshpts);
%-------------------------------------------------------%
% task-3. recover eigenpairs using SgEncd Green's function
%-------------------------------------------------------%
load('Data/SgEncdGreen_Fem_Value.mat');
% compute diagonal values via numerical integration
NN_G_diag_Qpts = reshape(NN_G_diag_Qpts, [N_y * N_y, num_pts])';
weights_diag = reshape(weights_diag, [N_y * N_y, num_pts])';
NN_G_diag = sum( weights_diag .* NN_G_diag_Qpts, 2) ./ (pi * eps ^ 2);
% generate predicted green's matrix on grids
NN_G_Fem = reshape(NN_G_Fem, num_pts, num_pts);
NN_G_Fem(1:size(NN_G_Fem,1)+1:end) = NN_G_diag;
% solve eigen problems
num_eigvalues = 300;
lambda_pred = eigensolver_green_poisson2d(node, elem+1, NN_G_Fem, num_eigvalues);
lambda_exact = generate_exact_eigenpairs_poisson2d(num_eigvalues);
plot_eigenpair_Green_Poisson2D(lambda_pred, lambda_exact, num_eigvalues);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











end