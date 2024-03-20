function fem = generateMassMatrix(fem, mesh, boundary_conditions)


%% assemble mass matrix with all elements
fem.Sd = fem.S .* reshape(fem.d, 3, 1, []);
fem.Sd = reshape(fem.Sd, 3, 1, fem.n_Gauss_points, []);
fem.Sct = conj(permute(fem.Sd, [2 1 3 4]));
fem.SS = fem.Sd .* fem.Sct;
fem.SS = reshape(fem.SS, 9, fem.n_Gauss_points, []);
fem.ele_mass = reshape(pagemtimes(fem.SS, fem.prod_tensor11_wxarea), 1, []).';
fem.mass = sparse(fem.row_indices, fem.col_indices, fem.ele_mass, mesh.n_nodes, mesh.n_nodes);

fem.Sd_0 = fem.S .* reshape(fem.d_0, 3, 1, []);
fem.Sd_0 = reshape(fem.Sd_0, 3, 1, fem.n_Gauss_points, []);
fem.Sct_0 = conj(permute(fem.Sd_0, [2 1 3 4]));
fem.SS_0 = fem.Sd_0 .* fem.Sct_0;
fem.SS_0 = reshape(fem.SS_0, 9, fem.n_Gauss_points, []);
fem.ele_mass_0 = reshape(pagemtimes(fem.SS_0, fem.prod_tensor11_wxarea), 1, []).';
fem.mass_0 = sparse(fem.row_indices, fem.col_indices, fem.ele_mass_0, mesh.n_nodes, mesh.n_nodes);

fem.Sd_end = fem.S .* reshape(fem.d_end, 3, 1, []);
fem.Sd_end = reshape(fem.Sd_end, 3, 1, fem.n_Gauss_points, []);
fem.Sct_end = conj(permute(fem.Sd_end, [2 1 3 4]));
fem.SS_end = fem.Sd_end .* fem.Sct_end;
fem.SS_end = reshape(fem.SS_end, 9, fem.n_Gauss_points, []);
fem.ele_mass_end = reshape(pagemtimes(fem.SS_end, fem.prod_tensor11_wxarea), 1, []).';
fem.mass_end = sparse(fem.row_indices, fem.col_indices, fem.ele_mass_end, mesh.n_nodes, mesh.n_nodes);

%% assemble mass matrix with all elements (mainly concern with the infinite opposite elements)
% if ismember('Dirichlet', boundary_conditions)
%     fem.Sd_infinite_opposite = fem.S .* reshape(fem.d_infinite_opposite, 3, 1, []);  % nothing to do with fem.d here??
%     fem.Sd_infinite_opposite = reshape(fem.Sd_infinite_opposite, 3, 1, fem.n_Gauss_points, []);
%     fem.Sct_infinite_opposite = conj(permute(fem.Sd_infinite_opposite, [2 1 3 4]));
%     fem.SS_infinite_opposite = fem.Sd_infinite_opposite .* fem.Sct_infinite_opposite;
%     fem.SS_infinite_opposite = reshape(fem.SS_infinite_opposite, 9, fem.n_Gauss_points, []);
%     fem.ele_mass_infinite_opposite = reshape(pagemtimes(fem.SS_infinite_opposite, fem.prod_tensor11_wxarea), 1, []).';
%     fem.mass_infinite_opposite = sparse(fem.row_indices, fem.col_indices, fem.ele_mass_infinite_opposite, mesh.n_nodes, mesh.n_nodes);
% 
%     fem.Sd_infinite = fem.S .* reshape(fem.d_infinite, 3, 1, []);  % nothing to do with fem.d here??
%     fem.Sd_infinite = reshape(fem.Sd_infinite, 3, 1, fem.n_Gauss_points, []);
%     fem.Sct_infinite = conj(permute(fem.Sd_infinite, [2 1 3 4]));
%     fem.SS_infinite = fem.Sd_infinite .* fem.Sct_infinite;
%     fem.SS_infinite = reshape(fem.SS_infinite, 9, fem.n_Gauss_points, []);
%     fem.ele_mass_infinite = reshape(pagemtimes(fem.SS_infinite, fem.prod_tensor11_wxarea), 1, []).';
%     fem.mass_infinite = sparse(fem.row_indices, fem.col_indices, fem.ele_mass_infinite, mesh.n_nodes, mesh.n_nodes);
% end

% fem.tensor11_purely_interior = kron(mesh.purely_interior_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + ... 
%     kron(ones(1, mesh.n_purely_interior_elements), 1);
% tensor11_purely_interior = reshape(fem.tensor11_purely_interior, 1, 1, mesh.n_purely_interior_elements);
% fem.area_purely_interior = fem.area(mesh.index.element_purely_interior);
% area_purely_interior = reshape(fem.area(mesh.index.element_purely_interior), 1, 1, mesh.n_purely_interior_elements);
% tmp = [2, 1, 1; 1, 2, 1; 1, 1, 2];
% ele_mass_purely_interior = repmat(tmp, 1, mesh.n_purely_interior_elements);
% ele_mass_purely_interior = reshape(ele_mass_purely_interior, 3, 3, mesh.n_purely_interior_elements);
% ele_mass_purely_interior = pagemtimes(ele_mass_purely_interior, ... 
%                                pagemtimes(area_purely_interior, tensor11_purely_interior));
% fem.ele_mass_purely_interior = reshape(ele_mass_purely_interior, 1, 3 * 3 * mesh.n_purely_interior_elements).';
% fem.mass_purely_interior = sparse(fem.purely_interior_row_indices, fem.purely_interior_col_indices, ... 
%     fem.ele_mass_purely_interior, mesh.n_nodes, mesh.n_nodes);
% 
% 
% % fem.tensor11_purely_interior = kron(mesh.purely_interior_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + ... 
% %     kron(ones(1, mesh.n_purely_interior_elements), 1);
% % r = fem.rspts(:, 1);
% % s = fem.rspts(:, 2);
% % x = mesh.nodes(1, :);
% % x = x(mesh.purely_interior_elements(1 : 3, :));
% % y = mesh.nodes(2, :);
% % y = y(mesh.purely_interior_elements(1 : 3, :));
% % eval(['[fem.S_purely_interior, fem.detJ_purely_interior] = isoparametricMapVectorize(x, y, r, s, @' parameters.fem.order ');']);
% % 
% % fem.wxarea_purely_interior = (fem.qwgts / 2) .* fem.detJ_purely_interior;
% % fem.SS_purely_interior = zeros(9, length(r));
% % for i = 1 : length(r)
% %     fem.SS_purely_interior(:, i) = reshape(fem.S_purely_interior(:, i) * fem.S_purely_interior(:, i)', 9, 1);
% % end
% % fem.ele_mass_purely_interior = reshape(fem.SS_purely_interior * (fem.tensor11_purely_interior .* fem.wxarea_purely_interior), 1, []).';
% % fem.mass_purely_interior = sparse(fem.purely_interior_row_indices, fem.purely_interior_col_indices, ... 
% %     fem.ele_mass_purely_interior, mesh.n_nodes, mesh.n_nodes);
% 
% 
% % %% assemble stiffness matrix with purely inf elements
% % fem.tensor11_pure_inf = kron(mesh.pure_inf_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + kron(ones(1, mesh.n_pure_inf_elements), 1);
% % r = fem.rspts(:, 1);
% % s = fem.rspts(:, 2);
% % x = mesh.nodes(1, :);
% % x = x(mesh.pure_inf_elements(1 : 3, :));
% % y = mesh.nodes(2, :);
% % y = y(mesh.pure_inf_elements(1 : 3, :));
% % % [fem.S, fem.detJ] = isoparametricMapVectorize(x, y, r, s, @P1);
% % eval(['[fem.S_pure_inf, fem.detJ_pure_inf] = isoparametricMapVectorize(x, y, r, s, @' parameters.fem.order ');']);
% % 
% % fem.wxarea_pure_inf = (fem.qwgts / 2) .* fem.detJ_pure_inf;
% % fem.SS_pure_inf = zeros(9, length(r));
% % for i = 1 : length(r)
% %     fem.SS_pure_inf(:, i) = reshape(fem.S_pure_inf(:, i) * fem.S_pure_inf(:, i)', 9, 1);
% % end
% % fem.ele_mass_pure_inf = reshape(fem.SS_pure_inf * (fem.tensor11_pure_inf .* fem.wxarea_pure_inf), 1, []).';
% % fem.mass_inf = sparse(fem.pure_inf_row_indices, fem.pure_inf_col_indices, ... 
% %     fem.ele_mass_pure_inf, mesh.n_nodes, mesh.n_nodes);
% % 
% % 
% % %% assemble stiffness matrix with purely inf2 elements
% % fem.tensor11_pure_inf2 = kron(mesh.pure_inf2_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + kron(ones(1, mesh.n_pure_inf2_elements), 1);
% % r = fem.rspts(:, 1);
% % s = fem.rspts(:, 2);
% % x = mesh.nodes(1, :);
% % x = x(mesh.pure_inf2_elements(1 : 3, :));
% % y = mesh.nodes(2, :);
% % y = y(mesh.pure_inf2_elements(1 : 3, :));
% % % [fem.S, fem.detJ] = isoparametricMapVectorize(x, y, r, s, @P1);
% % eval(['[fem.S_pure_inf2, fem.detJ_pure_inf2] = isoparametricMapVectorize(x, y, r, s, @' parameters.fem.order ');']);
% % 
% % fem.wxarea_pure_inf2 = (fem.qwgts / 2) .* fem.detJ_pure_inf2;
% % fem.SS_pure_inf2 = zeros(9, length(r));
% % for i = 1 : length(r)
% %     fem.SS_pure_inf2(:, i) = reshape(fem.S_pure_inf2(:, i) * fem.S_pure_inf2(:, i)', 9, 1);
% % end
% % fem.ele_mass_pure_inf2 = reshape(fem.SS_pure_inf2 * (fem.tensor11_pure_inf2 .* fem.wxarea_pure_inf2), 1, []).';
% % fem.mass_inf2 = sparse(fem.pure_inf2_row_indices, fem.pure_inf2_col_indices, ... 
% %     fem.ele_mass_pure_inf2, mesh.n_nodes, mesh.n_nodes);
% % 
% % fem.mass = fem.mass_pure_interior + fem.mass_inf + fem.mass_inf2;
% 
% 
% 
% % %% assemble mass matrix with nodes on inf boundary
% % fem.tensor11_inf = kron(mesh.inf_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + kron(ones(1, mesh.n_inf_elements), 1);
% % r = fem.rspts(:, 1);
% % s = fem.rspts(:, 2);
% % x = mesh.nodes(1, :);
% % x = x(mesh.inf_elements(1 : 3, :));
% % y = mesh.nodes(2, :);
% % y = y(mesh.inf_elements(1 : 3, :));
% % % [fem.S, fem.detJ] = isoparametricMapVectorize(x, y, r, s, @P1);
% % eval(['[fem.S_inf, fem.detJ_inf] = isoparametricMapVectorize(x, y, r, s, @' parameters.fem.order ');']);
% % 
% % fem.wxarea_inf = (fem.qwgts / 2) .* fem.detJ_inf;
% % fem.SS_inf = zeros(9, length(r));
% % for i = 1 : length(r)
% %     fem.SS_inf(:, i) = reshape(fem.S_inf(:, i) * fem.S_inf(:, i)', 9, 1);
% % end
% % fem.ele_mass_inf = reshape(fem.SS_inf * (fem.tensor11_inf .* fem.wxarea_inf), 1, []).';
% % fem.mass_inf = sparse(fem.inf_row_indices, fem.inf_col_indices, fem.ele_mass_inf, mesh.n_nodes, mesh.n_nodes);
% % 
% % % len_QBC2 = length(mesh.new.index_QBC2);
% % fem.mass_free_sort = fem.mass_free(mesh.new.index_free, mesh.new.index_free);
% % % fem.mass_free_small = fem.mass_free(1 : mesh.n_free_nodes, 1 : mesh.n_free_nodes);
% % fem.mass_inf_sort = fem.mass_inf(mesh.new.index, mesh.new.index);
% % fem.mass_inf_small = fem.mass_inf_sort(mesh.n_nodes - mesh.n_inf_nodes - mesh.n_near_inf_nodes + 1 : mesh.n_nodes, ... 
% %                                             mesh.n_nodes - mesh.n_inf_nodes - mesh.n_near_inf_nodes + 1 : mesh.n_nodes);
% % % figure, spy(fem.mass_free)
% % % figure, spy(fem.mass_inf)
% 
% 
% 
% % mass = sparse(mesh.n_nodes, mesh.n_nodes);
% % for i = 1 : mesh.n_elements
% %     ele_mass = zeros(3, 3);
% %     for j = 1 : length(fem.qwgts) % quadrature loop
% %         r = fem.rspts(j, 1); % r coordinate
% %         s = fem.rspts(j, 2); % s coordinate
% %         x = mesh.nodes(1, mesh.elements(1 : 3, i));
% %         y = mesh.nodes(2, mesh.elements(1 : 3, i));
% %         [S, dSdx, dSdy, detJ] = isoparametricMap(x', y', r, s, @P1); % map
% %         wxarea = fem.qwgts(j) * detJ / 2; % weight times det(J)
% %         tensor11 = parameters.math.tensor11_fun(mesh.elements(4, i));
% %         % norm((S * S') * (tensor11 * wxarea) - reshape(fem.SS(:, j) * fem.wxarea(j, i), 3, 3), 1)
% %         ele_mass = ele_mass + (S * S') * (tensor11 * wxarea); % compute and add integrand to MK
% %     end
% %     mass(mesh.elements(1 : 3, i), mesh.elements(1 : 3, i)) = mass(mesh.elements(1 : 3, i), mesh.elements(1 : 3, i)) + ele_mass;
% % end
% % norm(full(fem.mass - mass), 1)

end