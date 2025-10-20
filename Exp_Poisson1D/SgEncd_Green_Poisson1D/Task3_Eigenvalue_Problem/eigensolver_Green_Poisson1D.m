function [mu, phi] = eigensolver_Green_Poisson1D(F, basis_type, P_mesh, T_mesh, V_basis, T_basis, eig_num)

format short e

%% 1. FEM Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------%
% 1-1. generate information and mesh matrices
%-----------------------------------------------------------%
number_of_elements = size(T_mesh,2);    
number_of_nodes = size(V_basis,2);

if basis_type == 101 % linear elements
    number_of_basis = 2; 
elseif basis_type == 102 % quadratic elements
    number_of_basis = 3;  
end
%-----------------------------------------------------------%
% 1-2. quadrature's points and weights on refenrece interval 
%-----------------------------------------------------------%
number_of_quadrature_points = 8;
[quadrature_weight_reference_1D, quadrature_point_reference_1D] = Quadrature_Rule_Reference_1D(number_of_quadrature_points);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% (a) solve generalized eigenvalue problem of KLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a-1) assemble stiffness matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lBaG = zeros(number_of_basis,number_of_quadrature_points,number_of_elements);
H = zeros(number_of_basis,number_of_basis,number_of_elements);
K = zeros(number_of_nodes);
phi_qpts = zeros(number_of_elements*number_of_quadrature_points,eig_num);

% stiffness matrix
for n = 1 : number_of_elements

    % generate quadrature weights and points on local element
    vertices = P_mesh(:,T_mesh(:,n));        
    lower_bound = min(vertices(1),vertices(2));
    upper_bound = max(vertices(1),vertices(2));        
    [quadrature_weight_local_1D, quadrature_point_local_1D] = Quadrature_Rule_Local_1D(quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound);

    for i = 1 : number_of_basis
        for j = 1 : number_of_quadrature_points
            %lBaG(i,j,n) denotes the value corresponding to ith local base of nth element at local jth Gauss point  
            lBaG(i,j,n) = Local_Basis_1D(quadrature_point_local_1D(j), vertices, basis_type, i, 0);                     
        end
    end
        
    for k = 1 : number_of_basis
        H(k,:,n) = quadrature_weight_local_1D * (repmat(lBaG(k,:,n),number_of_basis,1).*lBaG(:,:,n))';
    end
    K(T_basis(:,n),T_basis(:,n)) = H(:,:,n) + K(T_basis(:,n),T_basis(:,n));   
    
end

M = K*F*K;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a-2) solve generalized eigenvalue problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi,D] = eig(M,K);
[mu,index] = sort(diag(D),'descend');
phi = phi(:,index);

% normalize eigenfunctions
for i = 1 : number_of_nodes
    Normalize_matrices(1,i) = norm(phi(:,i));
end
Normalize_matrices = repmat(Normalize_matrices,number_of_nodes,1);
phi = phi./Normalize_matrices;

for i = 1 : size(V_basis,2)
   phi(:,i) = phi(:,i) ./ max( abs(phi(:,i)) ); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1 : number_of_elements
    phi_qpts((n-1) * number_of_quadrature_points + 1 : n*number_of_quadrature_points,:) = lBaG(:,:,n)'*phi(T_basis(:,n),1:eig_num); 
end

end

