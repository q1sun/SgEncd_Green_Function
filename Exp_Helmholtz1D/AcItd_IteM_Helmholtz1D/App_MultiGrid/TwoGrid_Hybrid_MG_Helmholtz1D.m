function [eh_L2_HbMG, rh_L2_HbMG] = TwoGrid_Hybrid_MG_Helmholtz1D(Ah, fh, uh, R_list, NN_Preconditioner, maxIt, omega)

matrix_size=size(Ah,1);
h = 1 / (matrix_size+1);

% generate stiffness matrix on the coarset grid
R = R_list{end};
P = 2 * R';
for i = length(R_list) - 1 : -1 : 1
    R = R * R_list{i};
    P = (2 * transpose(R_list{i})) * P;
end

uh_HbMG = zeros(matrix_size, maxIt); % record iterative solution
eh_HbMG = zeros(matrix_size, maxIt); % record approximation error
eh_HbMG(:,1) = uh(:) - uh_HbMG(:,1);
rh_HbMG = zeros(matrix_size, maxIt); % record residual
rh_HbMG(:,1) = fh(:) - Ah * uh_HbMG(:,1);

for i = 1: maxIt - 1
    uh_HbMG(:,i+1) = TwoGrid_Hybrid_Multigrid_Vcycle(Ah, R, P, fh, NN_Preconditioner, uh_HbMG(:,i), omega);
    eh_HbMG(:,i+1) = uh(:) - uh_HbMG(:,i+1); 
    rh_HbMG(:,i+1) = fh(:) - Ah * uh_HbMG(:,i+1);
end

eh_L2_HbMG = zeros(1, maxIt);
rh_L2_HbMG = zeros(1, maxIt);

for l =  1 : maxIt      
    eh_L2_HbMG(:, l) = sum( eh_HbMG(:,l) .* eh_HbMG(:,l) .* h, 1) .^ (1/2);
    rh_L2_HbMG(:, l) = sum( rh_HbMG(:,l) .* rh_HbMG(:,l) .* h, 1) .^ (1/2);
end

% Multigrid V-clcye
function x = TwoGrid_Hybrid_Multigrid_Vcycle(Ah, R, P, fh, NN_Preconditioner, x0, omega)

% Ah        : The coefficient matrix on the finest level
% R         : The restriction operator from the finest level to the coarsest
% P         : The prolongation operator from the coarsest level to the finest
% fh        : The right hand side
	
	% Load coefficient matrix
	A0 = R * Ah * P;
    M = diag(diag(Ah));
    N = M - Ah; 
	
	% 2 Pre-smoothing steps
	x = x0 + h * NN_Preconditioner * (fh - Ah * x0);
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

	% Compute residual and transfer to coarse grid
	rh   = fh - Ah * x;
	rH = R * rh;
	
	% Solve coarse grid problem
	eH = A0 \ rH;
	
	% Transfer error to fine grid and correct
	x = x + P * eH;
	
	% 2 Post-smoothing steps
	x = x + h * NN_Preconditioner * (fh - Ah * x);
    x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

end

end

