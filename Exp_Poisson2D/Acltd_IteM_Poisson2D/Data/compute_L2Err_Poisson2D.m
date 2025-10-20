function [eh_L2,rh_L2]=compute_L2Err_Poisson2D(lBaG,P_mesh,T_mesh,index_intrr,Ah,fh,uh,Vh_item)

number_of_elements=size(T_mesh,2);
number_of_quadrature_points=7;

eh_L2 = zeros(1, size(Vh_item,2));
rh_L2 = zeros(1, size(Vh_item,2));

Uh = repmat(uh,1,size(eh_L2,2));
Eh_Hybrid = Uh - Vh_item;
Eh_Hybrid_plot = zeros(size(P_mesh,2), size(Vh_item,2));
Eh_Hybrid_plot(index_intrr,:) = Eh_Hybrid;

Rh_Hybrid_plot = zeros( size(P_mesh,2) , size(Vh_item,2));
Rh_Hybrid = repmat(fh,1,size(eh_L2,2)) - Ah * Vh_item;
Rh_Hybrid_plot(index_intrr,:) = Rh_Hybrid;

for n = 1 : number_of_elements
    
    % generate quadrature points for the n-th element
    vertices_triangle = P_mesh(1:2, T_mesh(1:3,n));
    [quadrature_weight_local_triangle,quadrature_point_local_triangle] = Quadrature_Rule_Local_Triangle(vertices_triangle', number_of_quadrature_points);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eh_L2 = eh_L2 + quadrature_weight_local_triangle * ((Eh_Hybrid_plot(T_mesh(:,n),:)'*lBaG(:,:,n))').^2;
    rh_L2 = rh_L2 + quadrature_weight_local_triangle * ((Rh_Hybrid_plot(T_mesh(:,n),:)'*lBaG(:,:,n))').^2;

end

eh_L2=sqrt(eh_L2);
rh_L2=sqrt(rh_L2);

end