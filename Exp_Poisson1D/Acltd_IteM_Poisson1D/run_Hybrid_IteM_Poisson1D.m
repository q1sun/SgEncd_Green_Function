function run_Hybrid_IteM_Poisson1D

format short e

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. subroutines
%-------------------------------------------------------%
addpath('Data\');
addpath('Plots\');
addpath('App_Hybrid_IteM');
%-------------------------------------------------------%
% 0-2. discrete system and save mesh points
%-------------------------------------------------------%
index_num = 8;
number_of_elements = 2 ^ index_num;
[Ah, fh, V_x_Gpts, eigVec_A] = discrete_system_Poisson1D(number_of_elements);
generate_meshPts_NN(V_x_Gpts, number_of_elements);
% calculate exact solution uh
uh = Ah \ fh;
%-------------------------------------------------------%
% 0-3. parameters of iteration
%-------------------------------------------------------%
omega = 0.5; % relaxation parameter
K = 2; L = 25;
number_of_iteration = K*L + 1;

%% 1. load neural preconditioner 
load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements));

%% 2. damped Jacobi method
[uh_Jacobi, eh_Jacobi, err_frq_Jacobi, ~, ~] = damped_Jacobi_Poisson1D(Ah, fh, uh, eigVec_A, number_of_iteration, omega);

%% 3. hybrid method using discrete neural Green's function
[uh_Hybrid, eh_Hybrid, err_frq_Hybrid, ~, ~] =hybrid_iteration_Poisson1D(Ah, fh, uh, eigVec_A, NN_Preconditioner, L, K, omega);

%% 4. plot comparitive results
% plot solution and error vector
    plot_Iter_U_Poisson1D(V_x_Gpts,uh_Hybrid,uh_Jacobi,uh,eh_Hybrid,eh_Jacobi,index_num,"");
% plot modewise error
    plot_FreErr_Poisson1D(err_frq_Hybrid',err_frq_Jacobi',index_num,""); % all modewise errors
    k=1;
    plot_FreErr_k_Poisson1D(err_frq_Hybrid,err_frq_Jacobi,k,index_num,"");
    k=floor(size(Ah,1)/2);
    plot_FreErr_k_Poisson1D(err_frq_Hybrid,err_frq_Jacobi,k,index_num,"");
    k=size(Ah,1);
    plot_FreErr_k_Poisson1D(err_frq_Hybrid,err_frq_Jacobi,k,index_num,"");
    
% different K values
[~, ~, ~, eh_L2_Jacobi, rh_L2_Jacobi] = damped_Jacobi_Poisson1D(Ah, fh, uh, eigVec_A, 201, omega);
[~, ~, ~, eh_L2_Hybrid_K2, rh_L2_Hybrid_K2] =hybrid_iteration_Poisson1D(Ah, fh, uh, eigVec_A, NN_Preconditioner, 100, 2, omega);
[~, ~, ~, eh_L2_Hybrid_K5, rh_L2_Hybrid_K5] =hybrid_iteration_Poisson1D(Ah, fh, uh, eigVec_A, NN_Preconditioner, 40, 5, omega);
[~, ~, ~, eh_L2_Hybrid_K10, rh_L2_Hybrid_K10] =hybrid_iteration_Poisson1D(Ah, fh, uh, eigVec_A, NN_Preconditioner, 20, 10, omega);
[~, ~, ~, eh_L2_Hybrid_K20, rh_L2_Hybrid_K20] =hybrid_iteration_Poisson1D(Ah, fh, uh, eigVec_A, NN_Preconditioner, 10, 20, omega);

plot_IterErr_Poisson1D(eh_L2_Hybrid_K2(:,1:201),rh_L2_Hybrid_K2(:,1:201),eh_L2_Jacobi(:,1:201),rh_L2_Jacobi(:,1:201),...
    eh_L2_Hybrid_K5(:,1:201), rh_L2_Hybrid_K5(:,1:201),...
eh_L2_Hybrid_K10(:,1:201), rh_L2_Hybrid_K10(:,1:201),...
eh_L2_Hybrid_K20(:,1:201), rh_L2_Hybrid_K20(:,1:201),...
index_num,"");




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end