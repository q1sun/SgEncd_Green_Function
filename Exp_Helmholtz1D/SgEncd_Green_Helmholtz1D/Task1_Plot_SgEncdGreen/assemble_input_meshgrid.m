function V_mesh = assemble_input_meshgrid(mesh_x, mesh_y)

% This function takes 1D coordinate vectors and creates a 2D mesh grid
% suitable for use as input to neural network models that require
% spatially distributed data.
%
% Inputs:
%   x_coords - 1D vector of x-coordinates
%   y_coords - 1D vector of y-coordinates
%
% Outputs:
%   mesh_grid - 2D grid structure containing:
%               .X - 2D matrix of x-coordinates
%               .Y - 2D matrix of y-coordinates
%               The grid spans all combinations of x_coords and y_coords

[mesh_x,mesh_y] = meshgrid(mesh_x,mesh_y);

V_mesh(1,:) = reshape(mesh_x,1,[]);
V_mesh(2,:) = reshape(mesh_y,1,[]);


end