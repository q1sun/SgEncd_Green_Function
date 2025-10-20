function [uh_Jacobi, eh_Jacobi, rh_Jacobi] = damped_Jacobi_Poisson2D(Ah, fh, uh, number_of_iteration, omega)

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split of matrix
M = diag(diag(Ah));
N = M - Ah; 
number_of_unknowns = size(Ah,1);

% damped Jacobi iteration
uh_Jacobi = zeros(number_of_unknowns, number_of_iteration); % record iterative solution

eh_Jacobi = zeros(number_of_unknowns, number_of_iteration); % record approximation error
eh_Jacobi(:,1) = uh(:) - uh_Jacobi(:,1);

rh_Jacobi = zeros(number_of_unknowns, number_of_iteration); % record residual
rh_Jacobi(:,1) = fh(:) - Ah * uh_Jacobi(:,1);

for i =  1 : number_of_iteration - 1
    uh_Jacobi(:,i+1) = omega * ( M^(-1) * N * uh_Jacobi(:,i) + M^(-1) * fh ) + (1-omega) * uh_Jacobi(:,i);
    eh_Jacobi(:,i+1) = uh(:) - uh_Jacobi(:,i+1); 
    rh_Jacobi(:,i+1) = fh(:) - Ah * uh_Jacobi(:,i+1);
    
end


end

