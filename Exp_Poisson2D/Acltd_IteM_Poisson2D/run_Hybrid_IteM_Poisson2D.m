function run_Hybrid_IteM_Poisson2D

format short e

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. subroutines
%-------------------------------------------------------%
addpath('Data\');
addpath('Plots\');
addpath('App_Hybrid_IteM\');
%-------------------------------------------------------%
% 0-2. discrete system and quadrature rule
%-------------------------------------------------------%
h1 = 0.2; % initial coarse mesh width for refinement process
index_num = 3; % no refinement
% generate the final mesh
[node, elem, isFreeNode] = generate_femgrids_disc2D(index_num, h1);
[Ah, fh, lBaG, ~] = discrete_system_Poisson2D(node', elem', isFreeNode');
uh = Ah \ fh;
%-------------------------------------------------------%
% 0-3. parameters of iteration
%-------------------------------------------------------%
omega = 0.5; % relaxation parameter
K = 2; L = 60;
number_of_iteration = K*L + 1;
%-------------------------------------------------------%
% 0-4. setup intput for pre-trained model
%-------------------------------------------------------%
freePts = node(isFreeNode, :);
weights_diag = prepare_NN_preconditioner(freePts, index_num);
%% 1. generate neural preconditioner 
load(sprintf('Data/NN_Preconditioner_h%d.mat',index_num));
NN_Preconditioner = generate_Neural_Preconditioner(freePts, weights_diag, NN_G_diag_Qpts, NN_G_Fem);

%% 2. damped Jacobi method
[uh_Jacobi, ~, ~] = damped_Jacobi_Poisson2D(Ah, fh, uh, number_of_iteration, omega);

%% 3. hybrid method using discrete neural Green's function
[uh_Hybrid, ~, ~] = hybrid_iteration_Poisson2D(Ah, fh, uh, NN_Preconditioner, L, K, omega);

%% 4. compute mode-wise errors
[err_frq_Hybrid, err_frq_Jacobi]=compute_FreErr_Poisson2D(Ah,omega,uh,uh_Hybrid,uh_Jacobi);
% plot modewise error
filename = sprintf("Meshsize_%d", index_num);
plot_FreErr_Poisson2D(err_frq_Hybrid,err_frq_Jacobi,filename,""); % all modewise errors
k=1;
plot_FreErr_k_Poisson2D(err_frq_Hybrid',err_frq_Jacobi',k,filename,"");
k=floor(size(Ah,1)/2);
plot_FreErr_k_Poisson2D(err_frq_Hybrid',err_frq_Jacobi',k,filename,"");
k=size(Ah,1);
plot_FreErr_k_Poisson2D(err_frq_Hybrid',err_frq_Jacobi',k,filename,"");

%% 5. compute L2-norm errors, residuals with various K values
[uh_Jacobi, ~, ~] = damped_Jacobi_Poisson2D(Ah, fh, uh, 201, omega);
[uh_Hybrid_K2, ~, ~] =hybrid_iteration_Poisson2D(Ah, fh, uh, NN_Preconditioner, 100, 2, omega);
[uh_Hybrid_K5, ~, ~] =hybrid_iteration_Poisson2D(Ah, fh, uh, NN_Preconditioner, 40, 5, omega);
[uh_Hybrid_K10, ~, ~] =hybrid_iteration_Poisson2D(Ah, fh, uh, NN_Preconditioner, 20, 10, omega);
[uh_Hybrid_K20, ~, ~] =hybrid_iteration_Poisson2D(Ah, fh, uh, NN_Preconditioner, 10, 20, omega);

[eh_L2_Jacobi, rh_L2_Jacobi]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_Jacobi);
[eh_L2_Hybrid_K2, rh_L2_Hybrid_K2]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_Hybrid_K2);
[eh_L2_Hybrid_K5, rh_L2_Hybrid_K5]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_Hybrid_K5);
[eh_L2_Hybrid_K10, rh_L2_Hybrid_K10]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_Hybrid_K10);
[eh_L2_Hybrid_K20, rh_L2_Hybrid_K20]=compute_L2Err_Poisson2D(lBaG,node',elem',isFreeNode',Ah,fh,uh,uh_Hybrid_K20);

plot_IterErr_Poisson2D(eh_L2_Hybrid_K2(:,1:201),rh_L2_Hybrid_K2(:,1:201),eh_L2_Jacobi(:,1:201),rh_L2_Jacobi(:,1:201),...
    eh_L2_Hybrid_K5(:,1:201), rh_L2_Hybrid_K5(:,1:201),...
eh_L2_Hybrid_K10(:,1:201), rh_L2_Hybrid_K10(:,1:201),...
eh_L2_Hybrid_K20(:,1:201), rh_L2_Hybrid_K20(:,1:201),...
filename,"");

end