function [eh_L2_MG, rh_L2_MG] = TwoGrid_MG_Poisson1D(Ah, fh, uh, R_list, maxIt, omega)

matrix_size=size(Ah,1);
h = 1 / (matrix_size+1);

% generate stiffness matrix on the coarset grid
R = R_list{end};
P = 2 * R';
for i = length(R_list) - 1 : -1 : 1
    R = R * R_list{i};
    P = (2 * transpose(R_list{i})) * P;
end

uh_MG = zeros(matrix_size, maxIt); % record iterative solution
eh_MG = zeros(matrix_size, maxIt); % record approximation error
eh_MG(:,1) = uh(:) - uh_MG(:,1);
rh_MG = zeros(matrix_size, maxIt); % record residual
rh_MG(:,1) = fh(:) - Ah * uh_MG(:,1);

for i = 1: maxIt - 1
    uh_MG(:,i+1) = TwoGrid_Multigrid_Vcycle(Ah, R, P, fh, uh_MG(:,i), omega);
    eh_MG(:,i+1) = uh(:) - uh_MG(:,i+1); 
    rh_MG(:,i+1) = fh(:) - Ah * uh_MG(:,i+1);
end

eh_L2_MG = zeros(1, maxIt);
rh_L2_MG = zeros(1, maxIt);

for l =  1 : maxIt      
    eh_L2_MG(:, l) = sum( eh_MG(:,l) .* eh_MG(:,l) .* h, 1) .^ (1/2);
    rh_L2_MG(:, l) = sum( rh_MG(:,l) .* rh_MG(:,l) .* h, 1) .^ (1/2);
end

% Multigrid V-clcye
function x = TwoGrid_Multigrid_Vcycle(Ah, R, P, fh, x0, omega)

% Ah        : The coefficient matrix on the finest level
% R         : The restriction operator from the finest level to the coarsest
% P         : The prolongation operator from the coarsest level to the finest
% fh        : The right hand side
% direct_n  : Threshold for directly solving A_list(level) * x = b
	
	% Load coefficient matrix
	A0 = R * Ah * P;
    M = diag(diag(Ah));
    N = M - Ah; 
	
	% 2 Pre-smoothing steps
	x = omega * ( M^(-1) * N * x0 + M^(-1) * fh ) + (1-omega) * x0;
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

	% Compute residual and transfer to coarse grid
	rh   = fh - Ah * x;
	rH = R * rh;
	
	% Solve coarse grid problem
	eH = A0 \ rH;
	
	% Transfer error to fine grid and correct
	x = x + P * eH;
	
	% 2 Post-smoothing steps
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;
    x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

end

end

