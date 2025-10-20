function run_Neural_Preconditioner_Poisson2D

format short e

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. subroutines
%-------------------------------------------------------%
addpath('Data\');
addpath('Plots\');
%-------------------------------------------------------%
% 0-2. discrete system and save mesh points
%-------------------------------------------------------%
h0 = 0.2; % initial coarse mesh width for refinement process
index_num = 2; % no refinement
[node, elem, isFreeNode] = generate_femgrids_disc2D(index_num, h0);
[Ah, fh, lBaG, ~] = discrete_system_Poisson2D(node', elem', isFreeNode');
uh = Ah \ fh;
%-------------------------------------------------------%
% 0-3. setup intput for pre-trained model
%-------------------------------------------------------%
freePts = node(isFreeNode, :);
weights_diag = prepare_NN_preconditioner(freePts, index_num);
%% 1. load neural Green's function for each task
%-------------------------------------------------------%
% task-1. acquire SgEncd Green's values on mesh points
%-------------------------------------------------------%
load(sprintf('Data/NN_Preconditioner_h%d.mat',index_num));
NN_Preconditioner = generate_Neural_Preconditioner(freePts, weights_diag, NN_G_diag_Qpts, NN_G_Fem);
%-------------------------------------------------------%
% task-2. draw eigenvalue spectra of BA
%-------------------------------------------------------%
plot_preconditioned_eigenvalues(Ah, NN_Preconditioner,index_num);
%-------------------------------------------------------%
% task-3. bicg solving linear algebraic system
%-------------------------------------------------------%
tol = 1e-16;
max_iter = 1e4;
%-------------------------------------------------------%
[x_bcg,fl1,rrl_bcg_h1,iter_bcg,rv1] = bicg(Ah,fh,tol,max_iter) ;
[x_prectd_bcg,fl2,rr1_prectd_bcg, iter_prectd_bcg,rv2] = bicg(Ah,fh,tol,max_iter,@(x,flag)NN_Preconditioner * x);

[err_bcg, res_bcg]=compute_L2Err_Poisson2D(lBaG, node', elem', isFreeNode', Ah, fh, uh, x_bcg);
[err_prectd_bcg, res_prectd_bcg]=compute_L2Err_Poisson2D(lBaG, node', elem', isFreeNode', Ah, fh, uh, x_prectd_bcg);
%-------------------------------------------------------%
fprintf('cg h=%e: Iter %d Err %e \n ', h0/2^(index_num-1), iter_bcg, err_bcg);
fprintf('prectd bcg h=%e: Iter %d Err %e \n', h0/2^(index_num-1), iter_prectd_bcg, err_prectd_bcg);



