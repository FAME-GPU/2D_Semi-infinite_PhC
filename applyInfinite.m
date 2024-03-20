function fem = applyInfinite(fem, mesh, parameters, bd)

%% assemble stiffness matrix with purely infinite elements

switch bd
    case 'Dirichlet'
        coeff = 0.0;
    case 'Neumann'
        coeff = 1.0;
    case 'Robin'
        coeff = 1.0;
end
fem.a_end(mesh.index.infinite_nodes_in_elements) = fem.a_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.b_end(mesh.index.infinite_nodes_in_elements) = fem.b_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.c_end(mesh.index.infinite_nodes_in_elements) = fem.c_end(mesh.index.infinite_nodes_in_elements) * coeff;
fem.d_end(mesh.index.infinite_nodes_in_elements) = fem.d_end(mesh.index.infinite_nodes_in_elements) * coeff;

% coeff = 1.0;
% fem.a(mesh.index.infinite_nodes_in_elements) = fem.a(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.b(mesh.index.infinite_nodes_in_elements) = fem.b(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.c(mesh.index.infinite_nodes_in_elements) = fem.c(mesh.index.infinite_nodes_in_elements) * coeff;
% fem.d(mesh.index.infinite_nodes_in_elements) = fem.d(mesh.index.infinite_nodes_in_elements) * coeff;

% fem.inv_tensor22_purely_infinite = kron(mesh.elements(4, mesh.index.element_purely_infinite) ~= 1, parameters.math.inv_tensor22 - eye(2)) + kron(ones(1, mesh.n_purely_infinite_elements), eye(2));
% fem.bc_purely_infinite = zeros(3, 2 * mesh.n_purely_infinite_elements);
% fem.bc_purely_infinite(:, 1 : 2 : 2 * mesh.n_purely_infinite_elements - 1) = fem.b(:, mesh.index.element_purely_infinite);
% fem.bc_purely_infinite(:, 2 : 2 : 2 * mesh.n_purely_infinite_elements) = fem.c(:, mesh.index.element_purely_infinite);
% 
% tmp_inv_tensor22_purely_infinite = reshape(fem.inv_tensor22_purely_infinite, 2, 2, mesh.n_purely_infinite_elements);
% tmp_bc = reshape(fem.bc_purely_infinite, 3, 2, mesh.n_purely_infinite_elements);
% tmp_new_bc = pagemtimes(tmp_inv_tensor22_purely_infinite, pagetranspose(tmp_bc));
% 
% ele_stiffness_purely_infinite = pagemtimes(conj(tmp_bc), tmp_new_bc);
% ele_stiffness_purely_infinite = pagemtimes(ele_stiffness_purely_infinite, reshape(fem.area(mesh.index.element_purely_infinite), 1, 1, mesh.n_purely_infinite_elements));
% fem.ele_stiffness_purely_infinite = reshape(ele_stiffness_purely_infinite, 1, 3 * 3 * mesh.n_purely_infinite_elements).';
% 
% purely_infinite_elements = mesh.purely_infinite_elements(1 : 3, :);
% fem.purely_infinite_row_indices = reshape(repmat(purely_infinite_elements, 3, 1), [], 1);
% fem.purely_infinite_col_indices = reshape(repmat(purely_infinite_elements(:), 1, 3)', 1, [])';
% fem.stiffness_purely_infinite = sparse(fem.purely_infinite_row_indices, fem.purely_infinite_col_indices, fem.ele_stiffness_purely_infinite, mesh.n_nodes, mesh.n_nodes);

end