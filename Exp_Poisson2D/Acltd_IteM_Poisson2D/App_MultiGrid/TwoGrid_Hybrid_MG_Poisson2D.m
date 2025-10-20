function uh_MG = TwoGrid_Hybrid_MG_Poisson2D(Ah, fh, R_list, P_list, NN_Preconditioner, maxIt, omega)

matrix_size = size(Ah, 1);
uh_MG = zeros(matrix_size, maxIt); % record iterative solution

R = R_list{end};
P = P_list{end-1};
for i = length(R_list)-1 : -1 : 3
    R = R_list{i} * R;
    P = P * P_list{i-1};
end

for i = 1: maxIt - 1
    uh_MG(:,i+1) = TwoGrid_Hybrid_Multigrid_Vcycle(Ah, R, P, fh, NN_Preconditioner, uh_MG(:,i), omega);
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
	x = x0 + NN_Preconditioner * (fh - Ah * x0);
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

	% Compute residual and transfer to coarse grid
	rh   = fh - Ah * x;
	rH = R * rh;
	
	% Solve coarse grid problem
	eH = A0 \ rH;
	
	% Transfer error to fine grid and correct
	x = x + P * eH;
	
	% 2 Post-smoothing steps
	x = x + NN_Preconditioner * (fh - Ah * x);
    x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

end

end



