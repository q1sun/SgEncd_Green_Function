function NN_Preconditioner = generate_Neural_Preconditioner(node, weights_diag, NN_G_diag_Qpts, NN_G_Fem)

num_pts = size(node,1);
N_y = 5; eps = 1e-2;

% compute diagonal values via numerical integration
NN_G_diag_Qpts = reshape(NN_G_diag_Qpts, [N_y * N_y, num_pts])';
weights_diag = reshape(weights_diag, [N_y * N_y, num_pts])';
NN_G_diag = sum( weights_diag .* NN_G_diag_Qpts, 2) ./ (pi * eps ^ 2);
% generate predicted green's matrix on grids
NN_Preconditioner = reshape(NN_G_Fem, num_pts, num_pts);
NN_Preconditioner(1:size(NN_Preconditioner,1)+1:end) = NN_G_diag';

end