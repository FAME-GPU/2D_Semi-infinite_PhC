function fem = generateStiffnessMatrix(fem, mesh, boundary_conditions)


%% assemble stiffness matrix with all elements
fem.new_bc = pagemtimes(fem.inv_tensor22, pagetranspose(fem.bc));
ele_stiffness = pagemtimes(conj(fem.bc), fem.new_bc);
ele_stiffness = pagemtimes(ele_stiffness, reshape(fem.area, 1, 1, mesh.n_elements));
fem.ele_stiffness = reshape(ele_stiffness, 1, 3 * 3 * mesh.n_elements).'; 
fem.stiffness = sparse(fem.row_indices, fem.col_indices, fem.ele_stiffness, mesh.n_nodes, mesh.n_nodes);

fem.new_bc_0 = pagemtimes(fem.inv_tensor22, pagetranspose(fem.bc_0));
ele_stiffness_0 = pagemtimes(conj(fem.bc_0), fem.new_bc_0);
ele_stiffness_0 = pagemtimes(ele_stiffness_0, reshape(fem.area, 1, 1, mesh.n_elements));
fem.ele_stiffness_0 = reshape(ele_stiffness_0, 1, 3 * 3 * mesh.n_elements).'; 
fem.stiffness_0 = sparse(fem.row_indices, fem.col_indices, fem.ele_stiffness_0, mesh.n_nodes, mesh.n_nodes);

fem.new_bc_end = pagemtimes(fem.inv_tensor22, pagetranspose(fem.bc_end));
ele_stiffness_end = pagemtimes(conj(fem.bc_end), fem.new_bc_end);
ele_stiffness_end = pagemtimes(ele_stiffness_end, reshape(fem.area, 1, 1, mesh.n_elements));
fem.ele_stiffness_end = reshape(ele_stiffness_end, 1, 3 * 3 * mesh.n_elements).'; 
fem.stiffness_end = sparse(fem.row_indices, fem.col_indices, fem.ele_stiffness_end, mesh.n_nodes, mesh.n_nodes);

%% assemble stiffness matrix with all elements (mainly concern with the infinite opposite elements)
% if ismember('Dirichlet', boundary_conditions)
%     fem.new_bc_infinite_opposite = pagemtimes(fem.inv_tensor22, pagetranspose(fem.bc_infinite_opposite));
%     ele_stiffness_infinite_opposite = pagemtimes(conj(fem.bc_infinite_opposite), fem.new_bc_infinite_opposite);
%     ele_stiffness_infinite_opposite = pagemtimes(ele_stiffness_infinite_opposite, reshape(fem.area, 1, 1, mesh.n_elements));
%     fem.ele_stiffness_infinite_opposite = reshape(ele_stiffness_infinite_opposite, 1, 3 * 3 * mesh.n_elements).';
%     fem.stiffness_infinite_opposite = sparse(fem.row_indices, fem.col_indices, fem.ele_stiffness_infinite_opposite, mesh.n_nodes, mesh.n_nodes);
% 
%     fem.new_bc_infinite = pagemtimes(fem.inv_tensor22, pagetranspose(fem.bc_infinite));
%     ele_stiffness_infinite = pagemtimes(conj(fem.bc_infinite), fem.new_bc_infinite);
%     ele_stiffness_infinite = pagemtimes(ele_stiffness_infinite, reshape(fem.area, 1, 1, mesh.n_elements));
%     fem.ele_stiffness_infinite = reshape(ele_stiffness_infinite, 1, 3 * 3 * mesh.n_elements).';
%     fem.stiffness_infinite = sparse(fem.row_indices, fem.col_indices, fem.ele_stiffness_infinite, mesh.n_nodes, mesh.n_nodes);
% end


% %% assemble stiffness matrix with purely inf elements
% fem.inv_tensor22_pure_inf = kron(mesh.pure_inf_elements(4, :) ~= 1, parameters.math.inv_tensor22 - eye(2)) + kron(ones(1, mesh.n_pure_inf_elements), eye(2));
% fem.pure_inf_bc = zeros(3, 2 * mesh.n_pure_inf_elements);
% fem.pure_inf_bc(:, 1 : 2 : 2 * mesh.n_pure_inf_elements - 1) = fem.b(:, mesh.index.element_pure_inf);
% fem.pure_inf_bc(:, 2 : 2 : 2 * mesh.n_pure_inf_elements) = fem.c(:, mesh.index.element_pure_inf);
% 
% tmp_inv_tensor22_pure_inf = reshape(fem.inv_tensor22_pure_inf, 2, 2, mesh.n_pure_inf_elements);
% tmp_bc = reshape(fem.pure_inf_bc, 3, 2, mesh.n_pure_inf_elements);
% tmp_new_bc = pagemtimes(tmp_inv_tensor22_pure_inf, pagetranspose(tmp_bc));
% % fem.new_bc = reshape(pagetranspose(tmp_new_bc), 3, 2 * mesh.n_elements);
% 
% ele_stiffness_pure_inf = pagemtimes(pagetranspose(tmp_new_bc), pagetranspose(tmp_bc));
% ele_stiffness_pure_inf = pagemtimes(ele_stiffness_pure_inf, reshape(fem.area(mesh.index.element_pure_inf), 1, 1, mesh.n_pure_inf_elements));
% fem.ele_stiffness_pure_inf = reshape(ele_stiffness_pure_inf, 1, 3 * 3 * mesh.n_pure_inf_elements).';
% 
% pure_inf_elements = mesh.pure_inf_elements(1 : 3, :);
% fem.pure_inf_row_indices = reshape(repmat(pure_inf_elements, 3, 1), [], 1);
% fem.pure_inf_col_indices = reshape(repmat(pure_inf_elements(:), 1, 3)', 1, [])';
% % stiffness = accumarray([row_indices, col_indices], ele_stiffness(:), [mesh.n_nodes, mesh.n_nodes], [], [], true);
% fem.stiffness_inf = sparse(fem.pure_inf_row_indices, fem.pure_inf_col_indices, fem.ele_stiffness_pure_inf, mesh.n_nodes, mesh.n_nodes);
% 
% %% assemble stiffness matrix with purely inf2 elements
% fem.inv_tensor22_pure_inf2 = kron(mesh.pure_inf2_elements(4, :) ~= 1, parameters.math.inv_tensor22 - eye(2)) + kron(ones(1, mesh.n_pure_inf2_elements), eye(2));
% fem.pure_inf_bc = zeros(3, 2 * mesh.n_pure_inf2_elements);
% fem.pure_inf_bc(:, 1 : 2 : 2 * mesh.n_pure_inf2_elements - 1) = fem.b(:, mesh.index.element_pure_inf2);
% fem.pure_inf_bc(:, 2 : 2 : 2 * mesh.n_pure_inf2_elements) = fem.c(:, mesh.index.element_pure_inf2);
% 
% tmp_inv_tensor22_pure_inf2 = reshape(fem.inv_tensor22_pure_inf, 2, 2, mesh.n_pure_inf2_elements);
% tmp_bc = reshape(fem.pure_inf_bc, 3, 2, mesh.n_pure_inf2_elements);
% tmp_new_bc = pagemtimes(tmp_inv_tensor22_pure_inf2, pagetranspose(tmp_bc));
% % fem.new_bc = reshape(pagetranspose(tmp_new_bc), 3, 2 * mesh.n_elements);
% 
% ele_stiffness_pure_inf2 = pagemtimes(pagetranspose(tmp_new_bc), pagetranspose(tmp_bc));
% ele_stiffness_pure_inf2 = pagemtimes(ele_stiffness_pure_inf2, reshape(fem.area(mesh.index.element_pure_inf2), 1, 1, mesh.n_pure_inf2_elements));
% fem.ele_stiffness_pure_inf2 = reshape(ele_stiffness_pure_inf2, 1, 3 * 3 * mesh.n_pure_inf2_elements).';
% 
% pure_inf2_elements = mesh.pure_inf2_elements(1 : 3, :);
% fem.pure_inf2_row_indices = reshape(repmat(pure_inf2_elements, 3, 1), [], 1);
% fem.pure_inf2_col_indices = reshape(repmat(pure_inf2_elements(:), 1, 3)', 1, [])';
% % stiffness = accumarray([row_indices, col_indices], ele_stiffness(:), [mesh.n_nodes, mesh.n_nodes], [], [], true);
% fem.stiffness_inf2 = sparse(fem.pure_inf2_row_indices, fem.pure_inf2_col_indices, fem.ele_stiffness_pure_inf2, mesh.n_nodes, mesh.n_nodes);

% fem.stiffness = fem.stiffness_pure_interior + fem.stiffness_inf + fem.stiffness_inf2;




% fem.stiffness_sort = fem.stiffness_pure_interior(mesh.new.index_free, mesh.new.index_free);
% % fem.stiffness_pure_interior_small = fem.stiffness_pure_interior(len_QBC2 + 1 : len_QBC2 + mesh.n_pure_interior_nodes, len_QBC2 + 1 : len_QBC2 + mesh.n_pure_interior_nodes);
% fem.stiffness_inf_sort = fem.stiffness_inf(mesh.new.index, mesh.new.index);
% fem.stiffness_inf_small = fem.stiffness_inf_sort(mesh.n_nodes - mesh.n_inf_nodes - mesh.n_near_pure_inf_nodes + 1 : mesh.n_nodes, ... 
%                                             mesh.n_nodes - mesh.n_pure_inf_nodes - mesh.n_near_pure_inf_nodes + 1 : mesh.n_nodes);
% figure, spy(fem.stiffness_pure_interior)
% figure, spy(fem.stiffness_pure_inf)



% stiffness2 = sparse(mesh.n_nodes, mesh.n_nodes);
% for i = 1 : mesh.n_elements
%     ele_stiffness = zeros(3, 3);
%     inv_tensor22 = parameters.math.inv_tensor22_fun(mesh.elements(4, i));
%     new_bc = inv_tensor22 * [fem.b(:, i), fem.c(:, i)].';
%     for ii = 1 : 3
%         for jj = 1 : 3
%             ele_stiffness(ii, jj) = new_bc(1, ii) * fem.b(jj, i) + new_bc(2, ii) * fem.c(jj, i);
%         end
%     end
%     ele_stiffness = ele_stiffness * fem.area(i);
%     stiffness2(mesh.elements(1 : 3, i), mesh.elements(1 : 3, i)) = stiffness2(mesh.elements(1 : 3, i), mesh.elements(1 : 3, i)) + ele_stiffness;
% end
% norm(full(stiffness2 - stiffness2), 1)

end