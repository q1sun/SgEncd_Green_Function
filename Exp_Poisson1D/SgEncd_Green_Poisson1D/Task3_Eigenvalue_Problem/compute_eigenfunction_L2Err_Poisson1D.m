function Phi_deep_L2Err = compute_eigenfunction_L2Err_Poisson1D(Phi_deep, number_of_eigenfunctions, P_mesh, T_mesh, T_basis, basis_type)

format short e

number_of_quadrature_points = 8;

Phi_deep_L2Err = zeros(1,number_of_eigenfunctions);
index_eigenfunc = linspace(1,number_of_eigenfunctions,number_of_eigenfunctions);

% eigenfunctions
phi_exact_func = @(x) sqrt(2)* sin( repmat(index_eigenfunc,number_of_quadrature_points,1).* pi .* x); % therefore, is also the relative L2 error

Phi_deep = Phi_deep * sqrt(2);
phi_deep = repmat(sign(Phi_deep(2,index_eigenfunc)),size(Phi_deep,1),1).* Phi_deep(:, index_eigenfunc);

% fem setting
number_of_elements = size(T_mesh,2); 

if basis_type == 101 % linear elements
    number_of_basis = 2; 
elseif basis_type == 102 % quadratic elements
    number_of_basis = 3;  
end

[quadrature_weight_reference_1D, quadrature_point_reference_1D] = Quadrature_Rule_Reference_1D(number_of_quadrature_points);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% (a) compute L2 and L_infty errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a-1) assemble stiffness matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lBaG = zeros(number_of_basis,number_of_quadrature_points,number_of_elements);

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
          
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a-2) compute error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : number_of_elements
    
    % generate quadrature weights and points on local element
    vertices = P_mesh(:,T_mesh(:,n));        
    lower_bound = min(vertices(1),vertices(2));
    upper_bound = max(vertices(1),vertices(2));        
    [quadrature_weight_local_1D, quadrature_point_local_1D] = Quadrature_Rule_Local_1D(quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Phi_exact=feval(phi_exact_func, repmat(quadrature_point_local_1D(:),1,number_of_all_eigenvalues));
    %sign_phi_exact=repmat(sign(Phi_exact(1,index_eigenfunc)),size(Phi_exact,1),1);
    err_qpt = feval(phi_exact_func, repmat(quadrature_point_local_1D(:),1,number_of_eigenfunctions)) ...
                                                    - (phi_deep(T_basis(:,n),:)'*lBaG(:,:,n))'; 
	Phi_deep_L2Err = Phi_deep_L2Err + quadrature_weight_local_1D * (err_qpt.^2);
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Phi_deep_L2Err = sqrt(Phi_deep_L2Err);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end