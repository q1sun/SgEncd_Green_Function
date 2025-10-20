function [uh_Jacobi, eh_Jacobi, err_frq_Jacobi, eh_L2_Jacobi, rh_L2_Jacobi] = damped_Jacobi_Poisson1D(Ah, fh, uh, eigVec_A, number_of_iteration, omega)

format short e

% Input:
% Ah, fh, uh, eigVec_A : Discrete system matrix, right-hand side, 
%                        exact solution, and eigenvectors of A
% NN_Preconditioner    : Pre-trained neural preconditioner
% omega                : Smoothing factor
% K, L                 : User-specified parameters

%% split of matrix
M = diag(diag(Ah));
N = M - Ah; 
matrix_size = size(eigVec_A,1);
h = 1 / (matrix_size+1);

% damped Jacobi iteration
uh_Jacobi = zeros(matrix_size, number_of_iteration); % record iterative solution

eh_Jacobi = zeros(matrix_size, number_of_iteration); % record approximation error
eh_Jacobi(:,1) = uh(:) - uh_Jacobi(:,1);

rh_Jacobi = zeros(matrix_size, number_of_iteration); % record residual
rh_Jacobi(:,1) = fh(:) - Ah * uh_Jacobi(:,1);

for i =  1 : number_of_iteration - 1
    uh_Jacobi(:,i+1) = omega * ( M^(-1) * N * uh_Jacobi(:,i) + M^(-1) * fh ) + (1-omega) * uh_Jacobi(:,i);
    eh_Jacobi(:,i+1) = uh(:) - uh_Jacobi(:,i+1); 
    rh_Jacobi(:,i+1) = fh(:) - Ah * uh_Jacobi(:,i+1);
    
end

%% log error per iteration
err_frq_Jacobi = zeros(matrix_size, number_of_iteration);
eh_L2_Jacobi = zeros(1, number_of_iteration);
rh_L2_Jacobi = zeros(1, number_of_iteration);

for l =  1 : number_of_iteration 
                
    err_frq_Jacobi(:, l) = abs(sum(repmat( eh_Jacobi(:,l), 1, matrix_size ) .*  eigVec_A .* h))'; 
    eh_L2_Jacobi(:, l) = sum( eh_Jacobi(:,l) .* eh_Jacobi(:,l) .* h, 1) .^ (1/2);
    rh_L2_Jacobi(:, l) = sum( rh_Jacobi(:,l) .* rh_Jacobi(:,l) .* h, 1) .^ (1/2) / norm(fh);
    
end


end