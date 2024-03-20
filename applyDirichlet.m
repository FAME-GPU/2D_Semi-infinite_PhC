function fem = applyDirichlet(fem, mesh, parameters)

%% apply Dirichlet boundary condition to the shape functions 
coeff = 0;

fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;

fem.a_end(mesh.index.infinite_nodes_in_elements) = fem.a_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.b_end(mesh.index.infinite_nodes_in_elements) = fem.b_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.c_end(mesh.index.infinite_nodes_in_elements) = fem.c_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.d_end(mesh.index.infinite_nodes_in_elements) = fem.d_end(mesh.index.infinite_nodes_in_elements) * coeff;

% bc_purely_infinite_oppo site = zeros(3, 2 * mesh.n_purely_infinite_opposite_elements);
% % important! determine shape functions with Dirichlet functions
% index_purely_infinite_opposite = ismember(mesh.purely_infinite_opposite_elements(1 : 3, :), mesh.index.node_purely_infinite_opposite);
% tmp_b_purely_infinite_opposite = fem.b(:, mesh.index.element_purely_infinite_opposite);
% tmp_b_purely_infinite_opposite(index_purely_infinite_opposite) = gamma_Dirichlet * tmp_b_purely_infinite_opposite(index_purely_infinite_opposite);
% fem.bc_purely_infinite_opposite(:, 1 : 2 : 2 * mesh.n_purely_infinite_opposite_elements - 1) = gamma_Dirichlet * tmp_b_purely_infinite_opposite;
% tmp_c_purely_infinite_opposite = fem.c(:, mesh.index.element_purely_infinite_opposite);
% tmp_c_purely_infinite_opposite(index_purely_infinite_opposite) = gamma_Dirichlet * tmp_c_purely_infinite_opposite(index_purely_infinite_opposite);
% fem.bc_purely_infinite_opposite(:, 2 : 2 : 2 * mesh.n_purely_infinite_opposite_elements) = gamma_Dirichlet * tmp_c_purely_infinite_opposite;
% fem.bc_purely_infinite_opposite(:, 1 : 2 : 2 * mesh.n_purely_infinite_opposite_elements - 1) = fem.b(:, mesh.index.element_purely_infinite_opposite);
% fem.bc_purely_infinite_opposite(:, 2 : 2 : 2 * mesh.n_purely_infinite_opposite_elements) = fem.c(:, mesh.index.element_purely_infinite_opposite);

end