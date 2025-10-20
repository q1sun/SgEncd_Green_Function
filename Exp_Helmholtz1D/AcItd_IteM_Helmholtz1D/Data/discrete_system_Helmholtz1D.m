function [Ah,fh, mesh, eigVec_B] = discrete_system_Helmholtz1D(number_of_elements)

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the Helmholtz's problem in one-dimension
%    - (c(x)u'(x))' - (k(x))^2 u(x) = f(x),  0 < x < 1
%         u(0) = u(1) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
left = 0; right = 1; % domain

% base functions
c_exact = @(x) (x-2) .^ 2; dc_exact = @(x) 2 * (x-2);
k_exact = @(x) 15 * sin(10 * x);
u_exact = @(x) 10 * x - 10 * x.^2 + 0.5 * sin( 20 * pi * x.^3 ); 
% first and second derivatives: du/dx, d^2u/d^2x
du_exact = @(x) 10 - 20 * x + 30 * pi * x.^2 .* cos(20 * pi * x.^3);
d2u_exact = @(x) -20 + 60 * pi * x .* cos(20 * pi * x.^ 3) - 1800 * pi^2 * x.^4 .* sin(20 * pi * x .^3);
% forcing term f
f_exact = @(x) - dc_exact(x) .* du_exact(x) - c_exact(x) .* d2u_exact(x) - k_exact(x) .^2 .* u_exact(x); % test example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. central difference scheme on equidistant mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 1-1. assemble coefficient matrix and load vector
%-------------------------------------------------------%
% generate mesh
h = 1/number_of_elements;
mesh = left : h : right;
mesh_mid = mesh + h/2;

% assemble matrix 
number_mesh_nodes = size(mesh,2); % including boundary values

D= [1/h^2 .* ( c_exact((mesh_mid(1:end-2))') + c_exact((mesh_mid(2:end-1))') ) - k_exact(mesh(2:end-1)').^2, ...
    (-1/h^2) .* c_exact((mesh_mid(1:end-2))') ,...
    (-1/h^2) .* c_exact((mesh_mid(2:end-1))')];

d=[0 1 -1];
Ah=spdiags(D,d,number_mesh_nodes-2,number_mesh_nodes-2);

% load right-hand-side vector
fh = f_exact(mesh)';

% equivalent linear system by removing boundary nodes
fh =fh( 2:number_mesh_nodes-1 );

% treat Dirichlet boundary conditions
fh(1,1) = fh(1,1) +  0;
fh(end,1) = fh(end,1) + 0;
%-------------------------------------------------------%
% 1-2. compute eigenvectors of damped Jacobi iterative matrix M (omega = 0.5)
%-------------------------------------------------------%
Ah = full(Ah);
D = diag(diag(Ah));
B = eye(size(Ah)) - 0.5 * inv(D) * Ah;
[eigVec_B, eigVal_B] = eig(B);
% sort and orthonormalize
[eigVal_B,index] = sort(abs(diag(eigVal_B)),'descend');
eigVec_B = eigVec_B(:,index);
[eigVec_B, ~] = gsog(eigVec_B); 


end