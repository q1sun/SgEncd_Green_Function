function [uh_MG] = FourGrid_MG_Poisson2D(Ah, fh, R_list, P_list, maxIt, N0, omega)

matrix_size = size(Ah, 1);
uh_MG = zeros(matrix_size, maxIt); % record iterative solution

for i = 1: maxIt - 1
    uh_MG(:,i+1) = Multigrid_Vcycle(1, Ah, R_list, P_list, fh, N0, uh_MG(:,i), omega);
end

% V-cycle
function x = Multigrid_Vcycle(level, Ah, R_list, P_list, fh, direct_n, x0, omega)

% level     : The current level (initial is 1)
% R_list    : The array of restriction operators on each level
% P_list    : The array of prolongation operators on each level
% fh         : The right hand side
% direct_n  : Threshold for directly solving A0 * x = b
	
	% Load coefficient matrix
    num_level = size(R_list, 1);
    M = diag(diag(Ah));
    N = M - Ah; 
	
	% If the problem is small enough, solve it directly
	n = size(fh, 1);
	if (n <= direct_n)
		x = Ah \ fh;
		return;
	end

	% 2 Pre-smoothing steps
	x = omega * ( M^(-1) * N * x0 + M^(-1) * fh ) + (1-omega) * x0;
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;
	
	% Load restriction operator and construct interpolation operator
	R = cell2mat(R_list(num_level + 1 - level, 1));
	P = cell2mat(P_list(num_level - level, 1));
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	r   = fh - Ah * x;
	r_H = R * r;
	
	% Solve coarse grid problem recursively
    x0  = zeros(coarse_n, 1);
    AH = R * Ah * P;
	e_H = Multigrid_Vcycle(level + 1, AH, R_list, P_list, r_H, direct_n, x0, omega);
	
	% Transfer error to fine grid and correct
	x = x + P * e_H;
	
	% 2 Post-smoothing steps
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;
    x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

end

end

