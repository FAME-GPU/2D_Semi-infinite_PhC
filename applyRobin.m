function fem = applyRobin(fem, mesh, parameters)

%% assemble stiffness matrix with purely Robin elements
coeff = parameters.math.boundary_condition.Robin;

fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.a_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.b_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.c_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;
fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) = fem.d_0(mesh.index.infinite_opposite_nodes_in_elements) * coeff;

fem.a_end(mesh.index.infinite_nodes_in_elements) = fem.a_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.b_end(mesh.index.infinite_nodes_in_elements) = fem.b_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.c_end(mesh.index.infinite_nodes_in_elements) = fem.c_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.d_end(mesh.index.infinite_nodes_in_elements) = fem.d_end(mesh.index.infinite_nodes_in_elements) * coeff;

% fem.inv_tensor22_purely_infinite_opposite = kron(mesh.elements(4, mesh.index.element_purely_infinite_opposite) ~= 1, parameters.math.inv_tensor22 - eye(2)) + kron(ones(1, mesh.n_purely_infinite_opposite_elements), eye(2));
% fem.bc_purely_infinite_opposite = zeros(3, 2 * mesh.n_purely_infinite_opposite_elements);
% fem.bc_purely_infinite_opposite(:, 1 : 2 : 2 * mesh.n_purely_infinite_opposite_elements - 1) = fem.b(:, mesh.index.element_purely_infinite_opposite);
% fem.bc_purely_infinite_opposite(:, 2 : 2 : 2 * mesh.n_purely_infinite_opposite_elements) = fem.c(:, mesh.index.element_purely_infinite_opposite);
% 
% tmp_inv_tensor22_purely_infinite_opposite = reshape(fem.inv_tensor22_purely_infinite_opposite, 2, 2, mesh.n_purely_infinite_opposite_elements);
% tmp_bc = reshape(fem.bc_purely_infinite_opposite, 3, 2, mesh.n_purely_infinite_opposite_elements);
% tmp_new_bc = pagemtimes(tmp_inv_tensor22_purely_infinite_opposite, pagetranspose(tmp_bc));
% 
% ele_stiffness_purely_infinite_opposite = pagemtimes(conj(tmp_bc), tmp_new_bc);
% ele_stiffness_purely_infinite_opposite = pagemtimes(ele_stiffness_purely_infinite_opposite, reshape(fem.area(mesh.index.element_purely_infinite_opposite), 1, 1, mesh.n_purely_infinite_opposite_elements));
% fem.ele_stiffness_purely_infinite_opposite = reshape(ele_stiffness_purely_infinite_opposite, 1, 3 * 3 * mesh.n_purely_infinite_opposite_elements).';
% 
% purely_infinite_opposite_elements = mesh.purely_infinite_opposite_elements(1 : 3, :);
% fem.purely_infinite_opposite_row_indices = reshape(repmat(purely_infinite_opposite_elements, 3, 1), [], 1);
% fem.purely_infinite_opposite_col_indices = reshape(repmat(purely_infinite_opposite_elements(:), 1, 3)', 1, [])';
% fem.stiffness_purely_infinite_opposite = sparse(fem.purely_infinite_opposite_row_indices, fem.purely_infinite_opposite_col_indices, fem.ele_stiffness_purely_infinite_opposite, mesh.n_nodes, mesh.n_nodes);

end