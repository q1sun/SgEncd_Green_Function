function V_mesh = assemble_input_meshgrid(mesh_x, mesh_y)

% This function takes 2D coordinate vectors and creates a 4D mesh grid
% suitable for use as input to neural network models that require
% spatially distributed data.
%
% Inputs:
%   mesh_x - 2D vector of x-coordinates
%   mesh_y - 2D vector of y-coordinates
%
% Outputs:
%   mesh_grid - 4D grid structure containing:
%               .X - 4D matrix of x-coordinates
%               .Y - 4D matrix of y-coordinates
%               The grid spans all combinations of x_coords and y_coords

mesh_x_repeat = repmat(mesh_x', size(mesh_y, 2), 1);
mesh_y_repeat = repelem(mesh_y', size(mesh_x, 2), 1);

V_mesh(1:2,:) = mesh_x_repeat';
V_mesh(3:4,:) = mesh_y_repeat';


end