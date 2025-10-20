function [V_mesh, V_Qpts, W_Qpts] = generate_mesh_quadrature_points(left, right, number_of_elements, number_of_quadrature_points)

format short e

% mesh 
V_mesh = linspace(left, right, number_of_elements+1);

% quadrature points
[quadrature_weight_reference_1D, quadrature_point_reference_1D] = Quadrature_Rule_Reference_1D(number_of_quadrature_points);

% quadrature points in each element
V_Qpts = zeros(number_of_elements * number_of_quadrature_points, 1); % coordinates
W_Qpts = zeros(number_of_elements * number_of_quadrature_points, 1); % weights

for n = 1 : number_of_elements

    % generate quadrature weights and points on local element       
    lower_bound = V_mesh(n);
    upper_bound = V_mesh(n+1);
    [quadrature_weight_local_1D, quadrature_point_local_1D] = Quadrature_Rule_Local_1D(quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound);

    V_Qpts(number_of_quadrature_points * (n-1) + 1 : number_of_quadrature_points * n, 1) = quadrature_point_local_1D';
    W_Qpts(number_of_quadrature_points * (n-1) + 1 : number_of_quadrature_points * n, 1) = quadrature_weight_local_1D';
    
end

end


function [quadrature_weight_local_1D, quadrature_point_local_1D] = Quadrature_Rule_Local_1D(quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound)

format short e

% generate (Gauss) quadrature weights and quadrature points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation
quadrature_weight_local_1D = (upper_bound-lower_bound)/2 * quadrature_weight_reference_1D;
quadrature_point_local_1D = (upper_bound-lower_bound)/2 * quadrature_point_reference_1D + (upper_bound+lower_bound)/2;

end


function [quadrature_weight_reference_1D, quadrature_point_reference_1D] = Quadrature_Rule_Reference_1D(number_of_quadrature_points)

format short e

% generate (Gauss) quadrature weights and quadrature points on the reference interval [-1,1]
if number_of_quadrature_points == 2
    quadrature_weight_reference_1D = [1,1];
    quadrature_point_reference_1D = [-1/sqrt(3),1/sqrt(3)];
elseif number_of_quadrature_points == 4
    quadrature_weight_reference_1D = [0.3478548451,0.3478548451,0.6521451549,0.6521451549];
    quadrature_point_reference_1D = [0.8611363116,-0.8611363116,0.3399810436,-0.3399810436];
elseif number_of_quadrature_points == 8
    quadrature_weight_reference_1D = [0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834];
    quadrature_point_reference_1D = [0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425];
end

end



