function [FreErr_Hybrid,FreErr_Jacobi]=compute_FreErr_Poisson2D(A,omega,uh_exact_int,Vh_Hybrid_int,Vh_Jacobi_int)

M=diag(diag(A));
G_J=diag(ones(1,size(A,1)))-omega*M^(-1)*A;

% eigenparis of discrete linear system
[eigF_G, eigV_G] = eig(G_J);

[eigV_G,index] = sort(diag(eigV_G),'descend');
eigF_G = eigF_G(:,index);
[eigF_G_orthonormal,R]=gsog(eigF_G);

Uh_exact=repmat(uh_exact_int,1,size(Vh_Hybrid_int,2));
Eh_Hybrid=Uh_exact-Vh_Hybrid_int;
Eh_Jacobi=Uh_exact-Vh_Jacobi_int;

FreErr_Hybrid=abs(Eh_Hybrid'*eigF_G_orthonormal)./sqrt(size(A,1)-1);
FreErr_Jacobi=abs(Eh_Jacobi'*eigF_G_orthonormal)./sqrt(size(A,1)-1);

end
