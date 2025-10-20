function [Ah,fh, mesh, eigVec_A] = discrete_system_Poisson1D(number_of_elements)

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the Poisson's problem in one-dimension
%    - u''(x) = f(x),  0 < x < 1
%         u(0) = u(1) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
left = 0; right = 1; % domain

% exact solutions for elliptic PDE
f_PDE_exact = @(x) - 60 * pi * x .* cos(20 * pi * x.^3) + 1800 * pi^2 * x.^4 .* sin( 20 * pi * x.^3) + 20; % test example 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. finite element scheme on equidistant mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 1-1. assemble stiffness matrix and load vector
%-------------------------------------------------------%
% generate mesh
h = 1/number_of_elements;
mesh = left : h : right;

% assemble matrix 
number_mesh_nodes = size(mesh,2); % including boundary values
D=[2*ones(number_mesh_nodes-2,1),-1*ones(number_mesh_nodes-2,1),-1*ones(number_mesh_nodes-2,1)];
d=[0 1 -1];
A=spdiags(D,d,number_mesh_nodes-2,number_mesh_nodes-2);
Ah = A/h;

% load right-hand-side vector
fh = f_PDE_exact(mesh)';

% equivalent linear system by removing boundary nodes
fh = h * fh(2:number_mesh_nodes-1);

%-------------------------------------------------------%
% 1-2. eigenpairs of stiffness matrix
%-------------------------------------------------------%
system_nodes = mesh( 2:number_mesh_nodes-1 );
number_system_nodes = size(system_nodes, 2);

% eigenparis of discrete linear system
Ah = full(Ah);
[eigVec_A, eigVal_A] = eig(Ah);

[eigVal_A,index] = sort(diag(eigVal_A),'ascend');
eigVec_A = eigVec_A(:,index);
eigVec_A = eigVec_A .* repmat( sign(eigVec_A(1,:)), number_system_nodes, 1 );
for i = 1 : number_system_nodes
   eigVec_A(:,i) = sqrt(2) * eigVec_A(:,i) ./ max( abs(eigVec_A(:,i)) ); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end