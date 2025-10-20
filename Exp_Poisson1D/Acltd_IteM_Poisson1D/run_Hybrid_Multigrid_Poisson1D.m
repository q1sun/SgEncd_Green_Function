function run_Hybrid_Multigrid_Poisson1D

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
index_num = 8;
number_of_elements = 2 ^ index_num;
[Ah, fh, V_x_Gpts, eigVec_A] = discrete_system_Poisson1D(number_of_elements);
generate_meshPts_NN(V_x_Gpts, number_of_elements);
uh = Ah \ fh;
%-------------------------------------------------------%
% 0-3. parameters of iteration
%-------------------------------------------------------%
omega = 0.5; % smoothing factor
number_of_iteration = 30;
index_final = 4;
direct_n = 2 ^ index_final;
[A_list, R_list, max_level] = Multigrid1D_Vcycle_GenMat(Ah, direct_n);

%% 1. load neural preconditioner 
load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements));

%% 2. 5-grid Multigrid Method
[eh_L2_5gMG, rh_L2_5gMG] = FiveGrid_MG_Poisson1D(A_list, fh, uh, R_list, number_of_iteration, direct_n, omega);

%% 3. 2-grid Multigrid Method
[eh_L2_2gMG, rh_L2_2gMG] = TwoGrid_MG_Poisson1D(Ah, fh, uh, R_list, number_of_iteration, omega);

%% 4. 2-grid Hybrid Multigrid Method
[eh_L2_HbMG, rh_L2_HbMG] = TwoGrid_Hybrid_MG_Poisson1D(Ah, fh, uh, R_list, NN_Preconditioner, number_of_iteration, omega);

%% 5. plot comparitive results
    plot_IterErr_MG_Poisson1D(eh_L2_5gMG, rh_L2_5gMG,...
        eh_L2_HbMG, rh_L2_HbMG, eh_L2_2gMG, rh_L2_2gMG,...
        index_num,"MG");

end