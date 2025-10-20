function weights_diag = prepare_NN_preconditioner(node, index_num)

%% generate diagonal samples for numerical integration
N_y = 5; % 5 gauss points 1d
eps = 1e-2; % radius of each ball
[diag_Qpts, weights_diag] = generate_diagonal_gauss_points(node, eps, N_y);

%% generate grid samples over the domain
mesh_x_repeat = repmat(node, size(node, 1), 1);
mesh_y_repeat = repelem(node, size(node, 1), 1);

mesh_Gpts(1:2,:) = mesh_x_repeat';
mesh_Gpts(3:4,:) = mesh_y_repeat';

save(sprintf('Data/mesh_Gpts_h%d.mat',index_num), 'diag_Qpts', 'mesh_Gpts');
end