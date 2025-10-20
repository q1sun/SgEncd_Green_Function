function [uh_Hybrid, eh_Hybrid, err_frq_Hybrid, eh_L2_Hybrid, rh_L2_Hybrid] =hybrid_iteration_Poisson1D(Ah,fh,uh,eigVec_A,NN_Preconditioner,K,L,omega)

format short e

% Input:
% Ah, fh, uh, eigVec_A : Discrete system matrix, right-hand side, 
%                        exact solution, and eigenvectors of A
% NN_Preconditioner    : Pre-trained neural preconditioner
% omega                : Smoothing factor
% K, L                 : User-specified parameters

% Output:
% uh_Hybrid     : Approximate solution at each iteration
% eh_Hybrid     : Error vector at each iteration
% err_frq_Hybrid: Error components across all frequencies
% eh_L2_Hybrid  : L2-norm of the error at each iteration
% rh_L2_Hybrid  : L2-norm of the residual at each iteration


%% split of matrix
matrix_size = size(eigVec_A,1);
h = 1 / (matrix_size+1);
number_of_iteration = K * L + 1;
M = diag(diag(Ah));
N = M - Ah; 

%% hybrid method iteration
uh_Hybrid=zeros(matrix_size,number_of_iteration);

eh_Hybrid = zeros(matrix_size, number_of_iteration); % record approximation error
eh_Hybrid(:,1) = uh(:) - uh_Hybrid(:,1);
rh_Hybrid = zeros(matrix_size, number_of_iteration); % record residual
rh_Hybrid(:,1) = fh(:) - Ah * uh_Hybrid(:,1);

for l=1:K*L
    if mod(l,L)==0
        % calcluate approximation error by neural preconditioner
        Rh = fh - Ah* uh_Hybrid(:,l);
        e_repNG_MeshPts = NN_Preconditioner* Rh;
        uh_Hybrid(:,l+1) = uh_Hybrid(:,l) + e_repNG_MeshPts;
    else
        % damped Jacobi iteration
        uh_Hybrid(:,l+1) = omega * ( M^(-1) * N * uh_Hybrid(:,l) + M^(-1) * fh ) + (1-omega) * uh_Hybrid(:,l);
    end

    eh_Hybrid(:,l+1) = uh(:) - uh_Hybrid(:,l+1); 
    rh_Hybrid(:,l+1) = fh(:) - Ah * uh_Hybrid(:,l+1);

end

%% log error per iteration

err_frq_Hybrid = zeros(matrix_size, number_of_iteration);
eh_L2_Hybrid = zeros(1, number_of_iteration);
rh_L2_Hybrid = zeros(1, number_of_iteration);

for l =  1 : number_of_iteration             
    err_frq_Hybrid(:, l) = abs(sum(repmat( eh_Hybrid(:,l), 1, matrix_size ) .*  eigVec_A .* h))'; 
    eh_L2_Hybrid(:, l) = sum( (eh_Hybrid(:,l)) .^2 .* h, 1) .^ (1/2);
    rh_L2_Hybrid(:, l) = sum( (rh_Hybrid(:,l)) .^2 .* h, 1) .^ (1/2) / norm(fh);
end

end


