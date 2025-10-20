function generate_meshPts_NN(x_Gpts, number_of_elements)

x_Gpts = x_Gpts(2:end-1);
mesh_Gpts = [ reshape(repmat( x_Gpts, size(x_Gpts,2), 1 ), [], 1 ), ...
                    repmat(x_Gpts', size(x_Gpts,2), 1 )];
mesh_Z = abs(mesh_Gpts(:,1) - mesh_Gpts(:,2));

mesh_Gpts = [mesh_Gpts, mesh_Z];
mesh_Gpts_symtry = [mesh_Gpts(:,2), mesh_Gpts(:,1), mesh_Z];

save(sprintf('Data/mesh_Gpts_h%d.mat',number_of_elements), 'mesh_Gpts', 'mesh_Gpts_symtry');

end