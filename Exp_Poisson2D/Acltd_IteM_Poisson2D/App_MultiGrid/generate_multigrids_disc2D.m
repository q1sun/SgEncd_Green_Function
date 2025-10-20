function [node, elem, isFreeNode, N0] = generate_multigrids_disc2D(num_level, h0)

addpath('Data\mesh_generation\');

% initial coarse level
[node,elem] = circlemesh(0,0,1,h0);
[fixedNode,bdEdge,isBdNode] = findboundary(elem);
isFreeNode = ~isBdNode;
N0 = size(node(isFreeNode,:), 1);

% refine the mesh
if num_level > 1
    for j = 1 : num_level - 1
        [node,elem] = uniformrefine(node,elem);
    end
end

[fixedNode,bdEdge,isBdNode] = findboundary(elem);
isFreeNode = ~isBdNode;

end