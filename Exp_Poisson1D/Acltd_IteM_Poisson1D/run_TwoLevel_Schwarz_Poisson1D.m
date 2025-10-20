function run_TwoLevel_Schwarz_Poisson1D

%% 0. Problem setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. subroutines
%-------------------------------------------------------%
addpath('Data\');
addpath('Plot\');
addpath('App_Preconditioner\');
%-------------------------------------------------------%
% 0-2. generate discrete system
%-------------------------------------------------------%
f_PDE_exact = @(x) - 60 * pi * x .* cos(20 * pi * x.^3) + 1800 * pi^2 * x.^4 .* sin( 20 * pi * x.^3) + 20; % test example
index = 16; N = 2 ^ index - 1; 
h = 1 / (N + 1); mesh = 0 : h : 1;
e = ones(N,1);
Ah = spdiags([-e 2*e -e], -1:1, N, N) / h;
fh = h * f_PDE_exact(mesh(2:end-1))';
%-------------------------------------------------------%
% 0-3. parameters of iteration
%-------------------------------------------------------%
num_sub = 40; 
overlap = 15; 

%% 1. Generate Local Matrix for Two-Level-Overlapping
sub_precs = cell(num_sub, 1);
idx_sets = cell(num_sub, 1);
sub_size = floor((N+1) / num_sub);

P_sets = cell(num_sub, 1);

for i = 1:num_sub
    left = max(1, (i - 1) * sub_size - overlap);
    right = min(N, i * sub_size + overlap);
    idx = (left:right)';
    idx_sets{i} = idx;
    Ri = sparse(1:length(idx), idx, 1, length(idx), N);
    A_loc = Ah(idx, idx);
    P_sets{i} = @(x) Ri' * (A_loc \(Ri *(Ah * x)));
    P_sets{i} = Ri' * (A_loc \(Ri *(Ah)));
    sub_precs{i} = @(r) A_loc \ r;
end

% additive preconditioner 
P_ad = sparse(size(P_sets{1},1), size(P_sets{1},2));
for i = 1:length(P_sets)
    P_ad = P_ad + P_sets{i};
end

%% 2. Coarse Level Neural Preconditioner
Nc = 2^(10) - 1;                              
load(sprintf('Data/NN_Preconditioner_h%d.mat',Nc+1));

% Construct R0 for 1D uniform grids (linear interpolation)
x_coarse = linspace(0,1,Nc+2); 
R0 = sparse(Nc, N);
x_fine = mesh(2:end-1);

for j = 2:Nc+2
    elem_finePts_idx = find(x_fine <= x_coarse(j) & x_fine > x_coarse(j-1));
    elem_finePts = x_fine(elem_finePts_idx);
    if j == 2
        R0(j-1,elem_finePts_idx) = (elem_finePts' - x_coarse(j-1))/(x_coarse(j) - x_coarse(j-1));
    elseif j >=3 && j <= Nc+1
        R0(j-1,elem_finePts_idx) = (elem_finePts' - x_coarse(j-1))/(x_coarse(j) - x_coarse(j-1));
        R0(j-2,elem_finePts_idx) = (x_coarse(j) - elem_finePts')/(x_coarse(j) - x_coarse(j-1));
    else
        R0(j-2,elem_finePts_idx) = (x_coarse(j) - elem_finePts')/(x_coarse(j) - x_coarse(j-1));
    end
end

%% 2. Compute eigenvalues of Additive Preconditioner
P_ad_NN = P_ad + R0' * (sparse(NN_Preconditioner) * (R0 *Ah));
P_ad_NN_func = @(x) P_ad_NN * x;

% calculate eigenvalues (largest and smallest)
opts.tol = 1e-6;
opts.maxit = 300;
opts.isreal = true;
opts.issym = false;

lambda_max = eigs(P_ad_NN_func, size(Ah,1), 1, 'lm', opts); 
lambda_min = eigs(P_ad_NN_func, size(Ah,1), 1, 'sm', opts); 

plot_additive_prectd_eigvalues(lambda_max, lambda_min, index);
kappa = lambda_max / lambda_min;
fprintf('number of elements：2^%d, kappa of BA = %.2e \n', index, kappa);

%% 3. Solve Linear Systems by Bicg
tol = 1e-10;
maxit = 50;

coarse_prec = @(r) R0' * (NN_Preconditioner * (R0 * r));   % M0^{-1} r
[u_prectd,flag,relres,iter] = bicg(Ah, fh, tol, maxit, @(x,flag)Additive_Preconditioning(idx_sets, sub_precs, num_sub, coarse_prec, x));
    
uh = Ah \ fh;
err_L2 = sum( (u_prectd - uh) .^2 .* h, 1) .^ (1/2);

fprintf('BICG done. Iter = %d, Err(L2) = %.2e \n', iter, err_L2);

end




