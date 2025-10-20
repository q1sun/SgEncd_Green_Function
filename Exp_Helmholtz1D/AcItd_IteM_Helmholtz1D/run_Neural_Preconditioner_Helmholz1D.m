function run_Neural_Preconditioner_Helmholz1D

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
number_of_elements_h6 = 2^6;
h6 = 1 / number_of_elements_h6;
[Ah_h6, fh_h6, V_x_Gpts_h6, eigVec_B] = discrete_system_Helmholtz1D(number_of_elements_h6);
generate_meshPts_NN(V_x_Gpts_h6, number_of_elements_h6);

number_of_elements_h8 = 2^8;
h8 = 1 / number_of_elements_h8;
[Ah_h8, fh_h8, V_x_Gpts_h8, eigVec_B] = discrete_system_Helmholtz1D(number_of_elements_h8); % here, eigenfunctions of differential operator instead of eigenvectors of A is used
generate_meshPts_NN(V_x_Gpts_h8, number_of_elements_h8);

number_of_elements_h10 = 2^10;
h10 = 1 / number_of_elements_h10;
[Ah_h10, fh_h10, V_x_Gpts_h10, eigVec_B] = discrete_system_Helmholtz1D(number_of_elements_h10);
generate_meshPts_NN(V_x_Gpts_h10, number_of_elements_h10);

number_of_elements_h12 = 2^12;
h12 = 1 / number_of_elements_h12;
[Ah_h12, fh_h12, V_x_Gpts_h12, eigVec_B] = discrete_system_Helmholtz1D(number_of_elements_h12);
generate_meshPts_NN(V_x_Gpts_h12, number_of_elements_h12);

%% 1. load neural Green's function for each task
%-------------------------------------------------------%
% task-1. acquire SgEncd Green's values on mesh points
%-------------------------------------------------------%
load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements_h6));
NN_Preconditioner_h6 = NN_Preconditioner;

load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements_h8));
NN_Preconditioner_h8 = NN_Preconditioner;

load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements_h10));
NN_Preconditioner_h10 = NN_Preconditioner;

load(sprintf('Data/NN_Preconditioner_h%d.mat',number_of_elements_h12));
NN_Preconditioner_h12 = NN_Preconditioner;
%-------------------------------------------------------%
% task-2. draw eigenvalue spectra of BA
%-------------------------------------------------------%
preconditioned_matrix_h6 = h6 * NN_Preconditioner_h6 * Ah_h6;
preconditioned_matrix_h8 = h8 * NN_Preconditioner_h8 * Ah_h8;
preconditioned_matrix_h10 = h10 * NN_Preconditioner_h10 * Ah_h10;
preconditioned_matrix_h12 = h12 * NN_Preconditioner_h12 * Ah_h12;

% plot eigenvalues
plot_preconditioned_eigvalues(preconditioned_matrix_h6, preconditioned_matrix_h8, preconditioned_matrix_h10, preconditioned_matrix_h12);
%-------------------------------------------------------%
% task-3. gmres solving linear algebraic system
%-------------------------------------------------------%
tol = 1e-16;
max_iter = 1e4;
%-------------------------------------------------------%
[x_gmres_h6,fl1,rr1_gmres_h6,iter_gmres_h6,rv1] = gmres(Ah_h6,fh_h6,[],tol,max_iter) ;
[x_prectd_gmres_h6,fl2,rr1_prectd_gmres_h6,iter_prectd_gmres_h6,rv2] = gmres(Ah_h6,fh_h6,[],tol,[],@(x)h6 * NN_Preconditioner_h6 * x);

uh_h6 = Ah_h6 \ fh_h6;
err_gmres_h6 = sum( (x_gmres_h6 - uh_h6) .^2 .* h6, 1) .^ (1/2);
err_prectd_gmres_h6 = sum( (x_prectd_gmres_h6 - uh_h6) .^2 .* h6, 1) .^ (1/2);
%-------------------------------------------------------%
[x_gmres_h8,fl1,rr1_gmres_h8,iter_gmres_h8,rv1] = gmres(Ah_h8,fh_h8,[],tol,max_iter) ;
[x_prectd_gmres_h8,fl2,rr1_prectd_gmres_h8,iter_prectd_gmres_h8,rv2] = gmres(Ah_h8,fh_h8,[],tol,[],@(x)h8 * NN_Preconditioner_h8 * x);

uh_h8 = Ah_h8 \ fh_h8;
err_gmres_h8 = sum( (x_gmres_h8 - uh_h8) .^2 .* h8, 1) .^ (1/2);
err_prectd_gmres_h8 = sum( (x_prectd_gmres_h8 - uh_h8) .^2 .* h8, 1) .^ (1/2);
%-------------------------------------------------------%
[x_gmres_h10,fl1,rr1_gmres_h10,iter_gmres_h10,rv1] = gmres(Ah_h10,fh_h10,[],tol,max_iter) ;
[x_prectd_gmres_h10,fl2,rr1_prectd_gmres_h10,iter_prectd_gmres_h10,rv2] = gmres(Ah_h10,fh_h10,[],tol,[],@(x)h10 * NN_Preconditioner_h10 * x);

uh_h10 = Ah_h10 \ fh_h10;
err_gmres_h10 = sum( (x_gmres_h10 - uh_h10) .^2 .* h10, 1) .^ (1/2);
err_prectd_gmres_h10 = sum( (x_prectd_gmres_h10 - uh_h10) .^2 .* h10, 1) .^ (1/2);
%-------------------------------------------------------%
[x_gmres_h12,fl1,rr1_gmres_h12,iter_gmres_h12,rv1] = gmres(Ah_h12,fh_h12,[],tol,max_iter) ;
[x_prectd_gmres_h12,fl2,rr1_prectd_gmres_h12,iter_prectd_gmres_h12,rv2] = gmres(Ah_h12,fh_h12,[],tol,[],@(x)h12 * NN_Preconditioner_h12 * x);

uh_h12 = Ah_h12 \ fh_h12;
err_gmres_h12 = sum( (x_gmres_h12 - uh_h12) .^2 .* h12, 1) .^ (1/2);
err_prectd_gmres_h12 = sum( (x_prectd_gmres_h12 - uh_h12) .^2 .* h12, 1) .^ (1/2);
%-------------------------------------------------------%
fprintf('gmres h=1/2^6: Iter %d Err %e \n ', iter_gmres_h6(2), err_gmres_h6);
fprintf('prectd gmres h=1/2^6: Iter %d Err %e \n',iter_prectd_gmres_h6(2), err_prectd_gmres_h6);
fprintf('gmres h=1/2^8: Iter %d Err %e \n', iter_gmres_h8(2), err_gmres_h8);
fprintf('prectd gmres h=1/2^8: Iter %d Err %e \n',iter_prectd_gmres_h8(2), err_prectd_gmres_h8);
fprintf('gmres h=1/2^10: Iter %d Err %e \n', iter_gmres_h10(2), err_gmres_h10);
fprintf('prectd gmres h=1/2^10: Iter %d Err %e \n',iter_prectd_gmres_h10(2), err_prectd_gmres_h10);
fprintf('gmres h=1/2^12: Iter %d Err %e \n', iter_gmres_h12(2), err_gmres_h12);
fprintf('prectd gmres h=1/2^12: Iter %d Err %e \n',iter_prectd_gmres_h12(2), err_prectd_gmres_h12);





