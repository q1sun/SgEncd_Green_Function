function [A, b, lBaG, eigVec_A] = discrete_system_Poisson2D(P_linear, T_linear, isFreeNode)

format short e

%% 0. Problem Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------%
% 0-1. 2D Poisson equation with Dirichlet boundary condition
%        - div( grad(u(x)) ) = f(x) in unit disc
%                  u(x) = 0 on bndry
%-----------------------------------------------------------%
a = @(x,y) 1 + 0 .* x .* y; % coefficient
f = @(x,y) 4 * exp( x.^2 + y.^2 - 1 ) .* (x.^2 + y.^2 + 1); % right-hand-side
%-----------------------------------------------------------%
% 0-2. tunable parameters
%-----------------------------------------------------------%
number_of_quadrature_points = 7;
basis_type = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Mesh Generation and Information Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh matrices are the same as that of linear elements
P_mesh = P_linear;
T_mesh = T_linear;
% finite elements
if basis_type == 1 % linear elements
    P_basis = P_linear;
    T_basis = T_linear;
    number_of_local_basis = 3;
elseif basis_type == 2 % quadratic elements
    P_basis = P_quadratic;
    T_basis = T_quadratic;
    number_of_local_basis = 6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 2. Assemble Stiffness Matrix and Load Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------%
% 2-1. assemble stiffness matrix
%-----------------------------------------------------------%
number_of_unknowns = size(P_basis,2);
number_of_elements = size(T_mesh,2);

lBaG = zeros(number_of_local_basis,number_of_quadrature_points,number_of_elements);
H = zeros(number_of_local_basis,number_of_local_basis,number_of_elements);
M = zeros(number_of_unknowns);

lBaG1 = zeros(number_of_local_basis,number_of_quadrature_points,number_of_elements);
H1 = zeros(number_of_local_basis,number_of_local_basis,number_of_elements);
A1 = zeros(number_of_unknowns);

lBaG2 = zeros(number_of_local_basis,number_of_quadrature_points,number_of_elements);
H2 = zeros(number_of_local_basis,number_of_local_basis,number_of_elements);
A2 = zeros(number_of_unknowns);

for n = 1 : number_of_elements % go through all elements
    
    % generate quadrature points for the n-th element
    vertices_triangle = P_mesh(1:2, T_mesh(1:3,n));
    [quadrature_weight_local_triangle,quadrature_point_local_triangle] = Quadrature_Rule_Local_Triangle(vertices_triangle', number_of_quadrature_points);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : number_of_local_basis
        for j = 1 : number_of_quadrature_points
            % lBaG(i,j,n) denotes the value corresponding to i-th local basis of n-th element at local j-th quadrature point  
            lBaG(i,j,n) = Basis_Local_Triangular2D(quadrature_point_local_triangle(j,1), quadrature_point_local_triangle(j,2), vertices_triangle, basis_type, i, 0, 0);                    
        end
    end

    for k = 1 : number_of_local_basis
        H(k,:,n) = quadrature_weight_local_triangle * (repmat(lBaG(k,:,n),number_of_local_basis,1).*lBaG(:,:,n))';
    end
    M(T_basis(:,n),T_basis(:,n)) = H(:,:,n) + M(T_basis(:,n),T_basis(:,n)); % \phi_i \phi_j  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : number_of_local_basis
        for j = 1 : number_of_quadrature_points
            % lBaG1(i,j,n) denotes the value corresponding to x-derivative of i-th local bais of n-th element at local j-th quadrature point  
            lBaG1(i,j,n) = Basis_Local_Triangular2D(quadrature_point_local_triangle(j,1), quadrature_point_local_triangle(j,2), vertices_triangle, basis_type, i, 1, 0);       
        end
    end
    
    for k = 1 : number_of_local_basis
        H1(k,:,n) = quadrature_weight_local_triangle * ((repmat(lBaG1(k,:,n),number_of_local_basis,1).*lBaG1(:,:,n))'.*repmat(a(quadrature_point_local_triangle(:,1),quadrature_point_local_triangle(:,2)),1,number_of_local_basis));
    end
    A1(T_basis(:,n),T_basis(:,n)) = H1(:,:,n) + A1(T_basis(:,n),T_basis(:,n)); % (grad_x \phi_i)*(grad_x \phi_j)             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : number_of_local_basis
        for j = 1 : number_of_quadrature_points
            % lBaG2(i,j,n) denotes the value corresponding to y-derivative of i-th local bais of n-th element at local j-th quadrature point   
            lBaG2(i,j,n) = Basis_Local_Triangular2D(quadrature_point_local_triangle(j,1), quadrature_point_local_triangle(j,2), vertices_triangle, basis_type, i, 0, 1);
        end
    end
    
    for k = 1 : number_of_local_basis
        H2(k,:,n) = quadrature_weight_local_triangle * ((repmat(lBaG2(k,:,n),number_of_local_basis,1).*lBaG2(:,:,n))'.*repmat(a(quadrature_point_local_triangle(:,1),quadrature_point_local_triangle(:,2)),1,number_of_local_basis));
    end
    A2(T_basis(:,n),T_basis(:,n)) = H2(:,:,n) + A2(T_basis(:,n),T_basis(:,n)); % (grad_y \phi_i)*(grad_y \phi_j)                     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
end

A = A1 + A2;
%-----------------------------------------------------------%
% 2-2. assemble load vector
%-----------------------------------------------------------%
b = M * f(P_basis(1,:),P_basis(2,:))';
%-----------------------------------------------------------%
% 2-3. treat Dirichlet boundary condition
%-----------------------------------------------------------%
% remove boundary points
A = A(isFreeNode, isFreeNode);
b = b(isFreeNode, :);

% eigenparis of discrete linear system
[eigVec_A, eigVal_A] = eig(A);
[eigVal_A,index] = sort(diag(eigVal_A),'ascend');
eigVec_A = eigVec_A(:,index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end