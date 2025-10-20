function [quadrature_weight_local_1D, quadrature_point_local_1D] = Quadrature_Rule_Local_1D(quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound)

format short e

% generate (Gauss) quadrature weights and quadrature points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation
quadrature_weight_local_1D = (upper_bound-lower_bound)/2 * quadrature_weight_reference_1D;
quadrature_point_local_1D = (upper_bound-lower_bound)/2 * quadrature_point_reference_1D + (upper_bound+lower_bound)/2;

end