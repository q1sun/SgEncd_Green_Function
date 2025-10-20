function [node, elem, isFreeNode] = generate_femgrids_disc2D(num_level, h0)

addpath('Data\mesh_generation\');

% initial coarse level
[node,elem] = circlemesh(0,0,1,h0);
% refine the mesh
if num_level > 1
    for j = 1 : num_level - 1
        [node,elem] = uniformrefine(node,elem);
    end
end

[fixedNode,bdEdge,isBdNode] = findboundary(elem);
isFreeNode = ~isBdNode;

end