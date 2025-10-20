function [R_list, P_list] = Multigrid2D_Vcycle_GenMat(elem, isFreeNode, N0)

% generate hierarchical mesh
[HB, NL, level] = HBstructure(elem,N0);
[P_list, R_list, freeNodej] = transferoperator(HB,NL,isFreeNode);

end