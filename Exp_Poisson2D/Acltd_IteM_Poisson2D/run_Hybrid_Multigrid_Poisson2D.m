function run_Hybrid_Multigrid_Poisson2D

format short e

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. subroutines
%-------------------------------------------------------%
addpath('Data\');
addpath('Plots\');
addpath('App_MultiGrid\');
%-------------------------------------------------------%
% 0-2. discrete system and quadrature rule
%-------------------------------------------------------%
h1 = 0.5; % initial coarse mesh width for refinement process
index_num = 4; % 4 counts of refinement
% generate the final mesh
[node, elem, isFreeNode, N0] = generate_multigrids_disc2D(index_num, h1);
[Ah, fh, lBaG, ~] = discrete_system_Poisson2D(node', elem', isFreeNode');
uh = Ah \ fh;
%-------------------------------------------------------%
% 0-3. parameters of iteration
%-------------------------------------------------------%
omega = 0.5; % smoothing factor
number_of_iteration = 31;
[R_list, P_list] = Multigrid2D_Vcycle_GenMat(elem, isFreeNode, N0);
%-------------------------------------------------------%
% 0-4. setup intput for pre-trained model
%-------------------------------------------------------%
freePts = node(isFreeNode, :);
weights_diag = prepare_2grid_NN_preconditioner(freePts, index_num);

%% 1. load neural preconditioner 
load(sprintf('App_MultiGrid/NN_Preconditioner_h%d.mat',index_num));
NN_Preconditioner = generate_2grid_Neural_Preconditioner(freePts, weights_diag, NN_G_diag_Qpts, NN_G_Fem);

%% 2. 4-grid Multigrid Method
uh_4gMG = FourGrid_MG_Poisson2D(Ah, fh, R_list, P_list, number_of_iteration, N0, omega);
[eh_L2_4gMG, rh_L2_4gMG]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_4gMG);

%% 3. 2-grid Multigrid Method
uh_2gMG = TwoGrid_MG_Poisson2D(Ah, fh, R_list, P_list, number_of_iteration, omega);
[eh_L2_2gMG, rh_L2_2gMG]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_2gMG);

%% 4. 2-grid Hybrid Multigrid Method
uh_HbMG = TwoGrid_Hybrid_MG_Poisson2D(Ah, fh, R_list, P_list, NN_Preconditioner, number_of_iteration, omega);
[eh_L2_HbMG, rh_L2_HbMG]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_HbMG);

%% 5. plot comparitive results
filename = sprintf("Meshsize_%d", index_num);
plot_IterErr_MG_Poisson2D(eh_L2_4gMG, rh_L2_4gMG,...
        eh_L2_HbMG, rh_L2_HbMG, eh_L2_2gMG, rh_L2_2gMG,...
        filename,"MG");

end

