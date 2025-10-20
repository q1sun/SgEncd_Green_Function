function [xy_pairs, weights] = generate_diagonal_gauss_points(x0_set, eps, N_y)
    % Generates Gauss integration points around diagonal points in 2D
    %
    % Inputs:
    %   x0_set - N_x × 2 matrix of center points [x0_1, x0_2]
    %   eps    - Radius of integration circles
    %   N_y    - Number of Gauss points per circle
    %
    % Outputs:
    %   xy_pairs - (N_x*N_y²) × 4 matrix of [x0, y] point pairs
    %   weights  - Corresponding integration weights

    N_x = size(x0_set, 1);
    
    % Generate angular and radial Gauss points
    [theta, weights_theta] = lgwt(N_y, 0, 2*pi);
    [r, weights_r] = lgwt(N_y, 0, eps);
    
    % Create 2D polar grid
    [R, Theta] = meshgrid(r, theta);
    [W_R, W_Theta] = meshgrid(weights_r, weights_theta);
    
    % Flatten to column vectors
    R = R(:); Theta = Theta(:);
    W_R = W_R(:); W_Theta = W_Theta(:);
    
    % Convert to Cartesian offsets
    y_offset_x = R .* cos(Theta);
    y_offset_y = R .* sin(Theta);
    weights_pt = R .* W_R .* W_Theta;
    
    % Replicate center points for all integration points
    x0_repeated = repelem(x0_set, N_y * N_y, 1);
    
    % Compute final coordinates and weights
    y_x = x0_repeated(:, 1) + repmat(y_offset_x, N_x, 1);
    y_y = x0_repeated(:, 2) + repmat(y_offset_y, N_x, 1);
    weights = repmat(weights_pt, N_x, 1);
    
    xy_pairs = [x0_repeated'; y_x'; y_y'];
end