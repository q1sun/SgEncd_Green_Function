function [uh_Hybrid, eh_Hybrid, rh_Hybrid] =hybrid_iteration_Poisson2D(Ah,fh,uh,NN_Preconditioner,L,K,omega)

format short e

%% split of matrix
M = diag(diag(Ah));
N = M - Ah; 
number_of_iteration=K*L+1;
number_of_unknowns=size(Ah,1);

%% hybrid method iteration
uh_Hybrid=zeros(number_of_unknowns,number_of_iteration);

eh_Hybrid = zeros(number_of_unknowns, number_of_iteration); % record approximation error
eh_Hybrid(:,1) = uh(:) - uh_Hybrid(:,1);
rh_Hybrid = zeros(number_of_unknowns, number_of_iteration); % record residual
rh_Hybrid(:,1) = fh(:) - Ah * uh_Hybrid(:,1);

for l=1:L*K

    if mod(l,K)==0
        % calcluate approximation error by Neural Preconditioner
        Rh = fh - Ah * uh_Hybrid(:,l);
        e_repNG_MeshPts = NN_Preconditioner * Rh;
        uh_Hybrid(:,l+1) = uh_Hybrid(:,l) + e_repNG_MeshPts;
    else
        % damped Jacobi iteration
        uh_Hybrid(:,l+1) = omega * ( M^(-1) * N * uh_Hybrid(:,l) + M^(-1) * fh ) + (1-omega) * uh_Hybrid(:,l);
    end

    eh_Hybrid(:,l+1) = uh(:) - uh_Hybrid(:,l+1); 
    rh_Hybrid(:,l+1) = fh(:) - Ah * uh_Hybrid(:,l+1);

end


end


