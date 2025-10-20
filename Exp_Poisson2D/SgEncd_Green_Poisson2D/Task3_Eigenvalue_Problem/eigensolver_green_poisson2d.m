function lambda = eigensolver_green_poisson2d(V_basis,T_basis,G,num_pairs)
% This function computes numerical eigenpairs for elliptic eigenvalue problems
% using a piecewise linear finite element method.
%
% Input:  V_basis - mesh grids
%         T_basis - vertex indices of each triangular element
%            G    - NN prediction on grids
%         num_pairs - number of required eigen values 
% Solver:  piecewise linear element 
% Output:  numerical eigenpairs
% Setting:  D = [0,1]*[0,1]
% reference: [2005 Schwab] FE for elliptic problems with sto coeff, Fig 6 (top left)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (0) problem setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Utils');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh: use comsol to generate mesh information matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear element
basis_type = 1;
number_of_basis = 3;
number_of_elements = size(T_basis,2);    
number_of_nodes = size(V_basis,2);
number_of_Gauss_pts = 7; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (1) solve generalized eigenvalue problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1-1) assemble stiffness matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lBaG = zeros(number_of_basis,number_of_Gauss_pts,number_of_elements);
H = zeros(number_of_basis,number_of_basis,number_of_elements);
K = zeros(number_of_nodes);

% stiffness matrix
for n = 1 : number_of_elements
    %generate Gauss points for n-th element
    vertices_triangle = V_basis(1:2,T_basis(1:3,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle] = generate_Gauss_point_local_triangle(number_of_Gauss_pts,vertices_triangle');
    
    for i = 1 : number_of_basis
        for j = 1 : number_of_Gauss_pts
            %lBaG(i,j,n) denotes the value corresponding to ith local base of nth element at local jth Gauss point  
            lBaG(i,j,n) = triangular_local_basis(Gauss_point_local_triangle(j,1),Gauss_point_local_triangle(j,2),vertices_triangle,basis_type,i,0,0);                     
        end
    end
        
    for k = 1 : number_of_basis
        H(k,:,n) = Gauss_coefficient_local_triangle * (repmat(lBaG(k,:,n),number_of_basis,1).*lBaG(:,:,n))';
    end
    K(T_basis(:,n),T_basis(:,n)) = H(:,:,n) + K(T_basis(:,n),T_basis(:,n));   
end

M = K*G*K;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1-2) solve generalized eigenvalue problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi,D] = eig(M,K);
[lambda,index] = sort(diag(D),'descend');
lambda = lambda(1:num_pairs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end