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
