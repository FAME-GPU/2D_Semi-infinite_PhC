function fem = applyNone(fem, mesh, parameters)

%% apply none boundary condition to the shape functions 
%  take Dirichlet boundary condition (0) as the end one

% fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
% fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
% fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
% fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;

% coeff = 0;
% 
% fem.a_end(mesh.index.infinite_nodes_in_elements) = fem.a_end(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.b_end(mesh.index.infinite_nodes_in_elements) = fem.b_end(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.c_end(mesh.index.infinite_nodes_in_elements) = fem.c_end(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.d_end(mesh.index.infinite_nodes_in_elements) = fem.d_end(mesh.index.infinite_nodes_in_elements) * coeff;

end