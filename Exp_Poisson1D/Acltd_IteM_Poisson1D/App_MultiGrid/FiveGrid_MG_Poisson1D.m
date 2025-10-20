function [eh_L2_MG, rh_L2_MG] = FiveGrid_MG_Poisson1D(A_list, fh, uh, R_list, maxIt, direct_n, omega)

Ah_fine = A_list{1};
matrix_size=size(Ah_fine,1);
h=1/(matrix_size+1);

uh_MG = zeros(matrix_size, maxIt); % record iterative solution
eh_MG = zeros(matrix_size, maxIt); % record approximation error
eh_MG(:,1) = uh(:) - uh_MG(:,1);
rh_MG = zeros(matrix_size, maxIt); % record residual
rh_MG(:,1) = fh(:) - Ah_fine * uh_MG(:,1);

for i = 1: maxIt - 1
    uh_MG(:,i+1) = Multigrid_Vcycle(1, A_list, R_list, fh, direct_n, uh_MG(:,i), omega);
    eh_MG(:,i+1) = uh(:) - uh_MG(:,i+1); 
    rh_MG(:,i+1) = fh(:) - Ah_fine * uh_MG(:,i+1);
end

eh_L2_MG = zeros(1, maxIt);
rh_L2_MG = zeros(1, maxIt);

for l =  1 : maxIt      
    eh_L2_MG(:, l) = sum( eh_MG(:,l) .* eh_MG(:,l) .* h, 1) .^ (1/2);
    rh_L2_MG(:, l) = sum( rh_MG(:,l) .* rh_MG(:,l) .* h, 1) .^ (1/2);
end

% V-cycle
function x = Multigrid_Vcycle(level, A_list, R_list, fh, direct_n, x0, omega)

% level     : The current level (initial is 1)
% A_list    : The array of coefficient matrices on each level
% R_list    : The array of restriction operators on each level
% fh         : The right hand side
% direct_n  : Threshold for directly solving A_list(level) * x = b
	
	% Load coefficient matrix
	Ah = cell2mat(A_list(level));
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
	R = cell2mat(R_list(level));
	P = 2 * R';
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	rh = fh - Ah * x;
	rH = R * rh;
	
	% Solve coarse grid problem recursively
    x0  = zeros(coarse_n, 1);
	eH = Multigrid_Vcycle(level + 1, A_list, R_list, rH, direct_n, x0, omega);
	
	% Transfer error to fine grid and correct
	x = x + P * eH;
	
	% 2 Post-smoothing steps
	x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;
    x = omega * ( M^(-1) * N * x + M^(-1) * fh ) + (1-omega) * x;

end

end

