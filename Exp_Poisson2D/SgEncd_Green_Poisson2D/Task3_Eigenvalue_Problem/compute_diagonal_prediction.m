function diagonal_values = compute_diagonal_prediction(f_values, weights, eps, N_x, N_y)

% reshape values of  (N_x, N_y * N_y) 形式
f_reshaped = reshape(f_values, [N_y * N_y, N_x])';
weights_reshaped = reshape(weights, [N_y * N_y, N_x])';

f_integral = sum( weights_reshaped .* f_reshaped, 2);
volumn_sphere = pi * eps ^ 2;

diagonal_values = f_integral ./ volumn_sphere;

end