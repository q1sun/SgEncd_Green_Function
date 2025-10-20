function lambda_exact = generate_exact_eigenpairs_poisson2d(num_pairs)
% This function computes exact eigenpairs for poisson equations in two
% dimensions.

%% (0) problem setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigen_list = [];
% Loop over m and n to collect enough eigenvalues
max_m = 100; max_n = 100;  % Try up to m, n = 100

%% (1) obtain eigen values of 2d Laplacian operator on a unit disc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 0:max_m
    jmn = besselzero(m, max_n, 1); % Get first max_n zeros of J_m
    for n = 1:max_n
        lambda = jmn(n)^2;
        if m == 0
            eigen_list = [eigen_list; m, n, lambda, jmn(n), 0];
        else
            eigen_list = [eigen_list; m, n, lambda, jmn(n), 1];
            eigen_list = [eigen_list; m, n, lambda, jmn(n), 2];
        end
    end
end

%% (3) obtain eigen values of the green's function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort by eigenvalue
eigen_list = sortrows(eigen_list, 3); 
eigen_list = eigen_list(1:num_pairs, :);
lambda_exact = eigen_list(:,3).^(-1);

end
