function fems = applyQuasiperiodic(fems, meshs, parameters, gamma)

geometries = fieldnames(meshs);
ngeom = length(geometries);

for gg = 1 : ngeom

    mesh = meshs.(geometries{gg});
    fem = fems.(geometries{gg});
    % boundary_conditions = parameters.math.(geometries{gg}).boundary_conditions;

    fem.a = fem.a_except_quasiperiodic;
    fem.b = fem.b_except_quasiperiodic;
    fem.c = fem.c_except_quasiperiodic;
    fem.d = fem.d_except_quasiperiodic;

    fem.a_0 = fem.a_except_quasiperiodic_0;
    fem.b_0 = fem.b_except_quasiperiodic_0;
    fem.c_0 = fem.c_except_quasiperiodic_0;
    fem.d_0 = fem.d_except_quasiperiodic_0;

    fem.a_end = fem.a_except_quasiperiodic_end;
    fem.b_end = fem.b_except_quasiperiodic_end;
    fem.c_end = fem.c_except_quasiperiodic_end;
    fem.d_end = fem.d_except_quasiperiodic_end;

    %% assemble stiffness matrix with quasiperiodic2 elements
    coeff = gamma;

    fem.a(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.b(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.c(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.d(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;

    fem.a_0(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_0(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.b_0(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_0(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.c_0(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_0(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.d_0(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_0(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;

    fem.a_end(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_end(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.b_end(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_end(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.c_end(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_end(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    fem.d_end(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_end(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;

    % if ismember('Dirichlet', boundary_conditions)
    %     fem.a_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.b_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.c_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.d_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    % 
    %     fem.a_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.b_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.c_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    %     fem.d_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
    % end

    fems.(geometries{gg}) = fem;

end

% if ismember('Dirichlet', parameters.math.boundary_condition.self)
%     fem.a_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.b_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.c_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.d_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_infinite_opposite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
% 
%     fem.a_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.a_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.b_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.b_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.c_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.c_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
%     fem.d_infinite(mesh.index.quasiperiodic2_nodes_in_elements) = fem.d_infinite(mesh.index.quasiperiodic2_nodes_in_elements) * coeff;
% end

% fem.inv_tensor22_quasiperiodic2 = kron(mesh.elements(4, mesh.index.element_quasiperiodic{2}) ~= 1, parameters.math.inv_tensor22 - eye(2)) + ... 
%                                   kron(ones(1, mesh.n_quasiperiodic_elements(2)), eye(2));
% fem.bc_quasiperiodic2 = zeros(3, 2 * mesh.n_quasiperiodic_elements(2));
% fem.bc_quasiperiodic2(:, 1 : 2 : 2 * mesh.n_quasiperiodic_elements(2) - 1) = fem.b(:, mesh.index.element_quasiperiodic{2});
% fem.bc_quasiperiodic2(:, 2 : 2 : 2 * mesh.n_quasiperiodic_elements(2)) = fem.c(:, mesh.index.element_quasiperiodic{2});
% 
% tmp_inv_tensor22_quasiperiodic2 = reshape(fem.inv_tensor22_quasiperiodic2, 2, 2, mesh.n_quasiperiodic_elements(2));
% tmp_bc2 = reshape(fem.bc_quasiperiodic2, 3, 2, mesh.n_quasiperiodic_elements(2));
% tmp_new_bc2 = pagemtimes(tmp_inv_tensor22_quasiperiodic2, pagetranspose(tmp_bc2));
% 
% ele_stiffness_quasiperiodic2 = pagemtimes(conj(tmp_bc2), tmp_new_bc2);
% ele_stiffness_quasiperiodic2 = pagemtimes(ele_stiffness_quasiperiodic2, reshape(fem.area(mesh.index.element_quasiperiodic), 1, 1, mesh.n_quasiperiodic_elements(2)));
% fem.ele_stiffness_quasiperiodic2 = reshape(ele_stiffness_quasiperiodic2, 1, 3 * 3 * mesh.n_quasiperiodic_elements(2)).';
% 
% quasiperiodic2_elements = mesh.quasiperiodic_elements{2}(1 : 3, :);
% fem.quasiperiodic2_row_indices = reshape(repmat(quasiperiodic2_elements, 3, 1), [], 1);
% fem.quasiperiodic2_col_indices = reshape(repmat(quasiperiodic2_elements(:), 1, 3)', 1, [])';
% fem.stiffness_quasiperiodic2 = sparse(fem.quasiperiodic2_row_indices, fem.quasiperiodic2_col_indices, fem.ele_stiffness_quasiperiodic2, mesh.n_nodes, mesh.n_nodes);
% 
% % fem.inv_tensor22_QBC2 = kron(mesh.QBC2_elements(4, :) ~= 1, parameters.math.inv_tensor22 - eye(2)) + ... 
% %     kron(ones(1, mesh.n_QBC2_elements), eye(2));
% % tmp_inv_tensor22_QBC2 = reshape(fem.inv_tensor22_QBC2, 2, 2, []);
% % fem.bc_QBC2 = zeros(3, 2 * mesh.n_QBC2_elements);
% % % important! determine shape functions with QBC
% % index_QBC2_hat = ismember(mesh.QBC2_elements(1 : 3, :), mesh.index.node_QBC2);
% % tmp_b_QBC2 = fem.b(:, mesh.index.element_QBC2);
% % tmp_b_QBC2(index_QBC2_hat) = gamma * tmp_b_QBC2(index_QBC2_hat);
% % fem.bc_QBC2(:, 1 : 2 : 2 * mesh.n_QBC2_elements - 1) = gamma * tmp_b_QBC2;
% % tmp_c_QBC2 = fem.c(:, mesh.index.element_QBC2);
% % tmp_c_QBC2(index_QBC2_hat) = gamma * tmp_c_QBC2(index_QBC2_hat);
% % fem.bc_QBC2(:, 2 : 2 : 2 * mesh.n_QBC2_elements) = gamma * tmp_c_QBC2;
% % tmp_bc_QBC2 = reshape(fem.bc_QBC2, 3, 2, mesh.n_QBC2_elements);
% % tmp_new_bc_QBC2 = pagemtimes(tmp_inv_tensor22_QBC2, pagetranspose(tmp_bc_QBC2));
% % 
% % ele_stiffness_QBC2 = pagemtimes(conj(tmp_bc_QBC2), tmp_new_bc_QBC2);
% % ele_stiffness_QBC2 = pagemtimes(ele_stiffness_QBC2, reshape(fem.area(mesh.index.element_QBC2), 1, 1, mesh.n_QBC2_elements));
% % fem.ele_stiffness_QBC2 = reshape(ele_stiffness_QBC2, 1, 3 * 3 * mesh.n_QBC2_elements).';
% % 
% % QBC2_elements = mesh.QBC2_elements(1 : 3, :);
% % fem.row_indices_QBC2 = reshape(repmat(QBC2_elements, 3, 1), [], 1);
% % fem.col_indices_QBC2 = reshape(repmat(QBC2_elements(:), 1, 3)', 1, [])';
% % fem.stiffness_QBC2 = sparse(fem.row_indices_QBC2, fem.col_indices_QBC2, fem.ele_stiffness_QBC2, mesh.n_nodes, mesh.n_nodes);
% 
% AA = fem.stiffness;
% AA(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) = ... 
%     AA(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) + ... 
%     fem.stiffness_QBC2(mesh.new.index_QBC2_plus_near_QBC2_all, mesh.new.index_QBC2_plus_near_QBC2_all);
% 
% mtx.A0 = AA(mesh.new.index_free, mesh.new.index_free);
% mtx.A = AA;
% mtx.A(mesh.new.index_inf2, mesh.new.index_inf2) = mtx.A(mesh.new.index_inf2, mesh.new.index_inf2) + ... 
%     AA(mesh.new.index_inf, mesh.new.index_inf);
% mtx.A = mtx.A(mesh.new.index_free, mesh.new.index_free);
% mtx.Ae = AA; 
% mtx.Ae(mesh.new.index_inf2, mesh.new.index_inf2) = mtx.Ae(mesh.new.index_inf2, mesh.new.index_inf2) + ... 
%     AA(mesh.new.index_inf, mesh.new.index_inf);
% mtx.Ae = mtx.Ae(mesh.new.index_free_end, mesh.new.index_free_end);
% 
% EE = sparse(mesh.n_nodes, mesh.n_nodes);
% EE(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) = ... 
%     AA(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) - ... 
%     blkdiag(AA(mesh.new.new_index_near_inf, mesh.new.new_index_near_inf), ... 
%                         AA(mesh.new.new_index_inf, mesh.new.new_index_inf));
% mtx.E = sparse(mesh.n_free_nodes, mesh.n_free_nodes);
% n_row_E = length(mesh.new.index_inf);
% n_col_E = length(mesh.new.index_near_inf);
% mtx.Eh = EE(mesh.new.index_inf, mesh.new.index_near_inf);
% mtx.E(1 : n_row_E, end - n_col_E + 1 : end) = mtx.Eh;
% mtx.Ee = sparse(mesh.n_end_free_nodes, mesh.n_end_free_nodes);
% mtx.Ee(1 : n_row_E, end - n_col_E + 1 : end) = mtx.Eh;
% 
% % EE = fem.stiffness_inf;
% % EE(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) = ... 
% %     EE(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) + ... 
% %     fem.stiffness_QBC2(mesh.new.index_QBC2_plus_near_QBC2_all, mesh.new.index_QBC2_plus_near_QBC2_all);
% % EE(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) = ... 
% %     AA(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) - ... 
% %     blkdiag(AA(mesh.new.new_index_near_inf, mesh.new.new_index_near_inf), ... 
% %                         AA(mesh.new.new_index_inf, mesh.new.new_index_inf));
% % mtx.EE = sparse(mesh.n_free_nodes, mesh.n_free_nodes);
% % n_row_EE = length(mesh.new.index_inf);
% % n_col_EE = length(mesh.new.index_near_inf);
% % mtx.EEh = EE(mesh.new.index_inf, mesh.new.index_near_inf);
% % mtx.EE(1 : n_row_EE, end - n_col_EE + 1 : end) = mtx.EEh;
% 
% % mtx.E = EE(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf);
% % mtx.E = mtx.E - blkdiag(AA(mesh.new.new_index_near_inf, mesh.new.new_index_near_inf), ... 
% %                         AA(mesh.new.new_index_inf, mesh.new.new_index_inf));
% % mtx.E = mtx.E(length(mesh.new.new_index_inf) + 1 : end, 1 : length(mesh.new.new_index_near_inf));
% % mtx.E = 1;
% 
% 
% % mtx.A(mesh.new.index_QBC_plus_near_QBC2, mesh.new.index_QBC_plus_near_QBC2) = ... 
% %     mtx.A(mesh.new.index_QBC_plus_near_QBC2, mesh.new.index_QBC_plus_near_QBC2) + ... 
% %     fem.stiffness_QBC2(mesh.new.index_QBC2_plus_near_QBC2, mesh.new.index_QBC2_plus_near_QBC2);
% % mtx.A = mtx.A(mesh.new.index_free, mesh.new.index_free);
% % 
% % mtx.E = mtx.E;
% 
% 
% %% apply QBC to mass matrix
% fem.tensor11_QBC2 = kron(mesh.QBC2_elements(4, :) ~= 1, parameters.math.tensor11 - 1) + ... 
%     kron(ones(1, mesh.n_QBC2_elements), 1);
% r = fem.rspts(:, 1);
% s = fem.rspts(:, 2);
% x = mesh.nodes(1, :);
% x = x(mesh.QBC2_elements(1 : 3, :));
% y = mesh.nodes(2, :);
% y = y(mesh.QBC2_elements(1 : 3, :));
% % [fem.S, fem.detJ] = isoparametricMapVectorize(x, y, r, s, @P1);
% eval(['[fem.S_QBC2, fem.detJ_QBC2] = isoparametricMapVectorize(x, y, r, s, @' parameters.fem.order ');']);
% 
% fem.wxarea_QBC2 = (fem.qwgts / 2) .* fem.detJ_QBC2;
% fem.SS_QBC2 = zeros(9, length(r));
% for i = 1 : length(r)
%     fem.SS_QBC2(:, i) = reshape(fem.S_QBC2(:, i) * fem.S_QBC2(:, i)', 9, 1);
% end
% fem.ele_mass_QBC2 = reshape(fem.SS_QBC2 * (fem.tensor11_QBC2 .* fem.wxarea_QBC2), 1, []).';
% fem.mass_QBC2 = sparse(fem.row_indices_QBC2, fem.col_indices_QBC2, fem.ele_mass_QBC2, mesh.n_nodes, mesh.n_nodes);
% 
% BB = fem.mass;
% BB(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) = ... 
%     BB(mesh.new.index_QBC_plus_near_QBC2_all, mesh.new.index_QBC_plus_near_QBC2_all) + ... 
%     fem.mass_QBC2(mesh.new.index_QBC2_plus_near_QBC2_all, mesh.new.index_QBC2_plus_near_QBC2_all);
% 
% mtx.B0 = BB(mesh.new.index_free, mesh.new.index_free);
% mtx.B = BB;
% mtx.B(mesh.new.index_inf2, mesh.new.index_inf2) = mtx.B(mesh.new.index_inf2, mesh.new.index_inf2) + ... 
%     BB(mesh.new.index_inf, mesh.new.index_inf);
% mtx.B = mtx.B(mesh.new.index_free, mesh.new.index_free);
% mtx.Be = BB;
% mtx.Be(mesh.new.index_inf2, mesh.new.index_inf2) = mtx.Be(mesh.new.index_inf2, mesh.new.index_inf2) + ... 
%     BB(mesh.new.index_inf, mesh.new.index_inf);
% mtx.Be = mtx.Be(mesh.new.index_free_end, mesh.new.index_free_end);
% 
% 
% FF = sparse(mesh.n_nodes, mesh.n_nodes);
% FF(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) = ... 
%     BB(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf) - ... 
%     blkdiag(BB(mesh.new.new_index_near_inf, mesh.new.new_index_near_inf), ... 
%                         BB(mesh.new.new_index_inf, mesh.new.new_index_inf));
% mtx.F = sparse(mesh.n_free_nodes, mesh.n_free_nodes);
% n_row_F = length(mesh.new.index_inf);
% n_col_F = length(mesh.new.index_near_inf);
% mtx.Fh = FF(mesh.new.index_inf, mesh.new.index_near_inf);
% mtx.F(1 : n_row_F, end - n_col_F + 1 : end) = mtx.Fh;
% mtx.Fe = sparse(mesh.n_end_free_nodes, mesh.n_end_free_nodes);
% mtx.Fe(1 : n_row_E, end - n_col_E + 1 : end) = mtx.Eh;
% 
% % mtx.F = BB(mesh.new.index_near_inf_plus_inf, mesh.new.index_near_inf_plus_inf);
% % mtx.F = mtx.F - blkdiag(BB(mesh.new.new_index_near_inf, mesh.new.new_index_near_inf), ... 
% %                         BB(mesh.new.new_index_inf, mesh.new.new_index_inf));
% % mtx.F = mtx.F(length(mesh.new.new_index_inf) + 1 : end, 1 : length(mesh.new.new_index_near_inf));
% % mtx.F = 1;
% 
% 
% 
% % mtx.B(mesh.new.index_QBC_plus_near_QBC2, mesh.new.index_QBC_plus_near_QBC2) = ... 
% %     mtx.B(mesh.new.index_QBC_plus_near_QBC2, mesh.new.index_QBC_plus_near_QBC2) + ... 
% %     fem.mass_QBC2(mesh.new.index_QBC2_plus_near_QBC2, mesh.new.index_QBC2_plus_near_QBC2);
% % mtx.B = mtx.B(mesh.new.index_free, mesh.new.index_free);
% 
% % mtx.B = mtx.B + fem.mass_QBC2_small;
% 
% 
% % for i = 1 : size(mesh.QBC2_elements, 2)
% % 
% %     QBC_i = mesh.index.element_QBC(i);
% % 
% %     % add QBC for stiffness matrix
% %     ele_stiffness_QBC = zeros(3, 3);
% %     inv_tensor22_QBC = parameters.math.inv_tensor22_fun(mesh.QBC2_elements(4, i));
% %     new_bc = inv_tensor22_QBC * [b(:, QBC_i), c(:, QBC_i)].';
% %     for ii = 1 : 3
% %         for jj = 1 : 3
% %             ele_stiffness_QBC(ii, jj) = new_bc(1, ii) * b(jj, QBC_i) + new_bc(2, ii) * c(jj, QBC_i);
% %         end
% %     end
% %     ele_stiffness_QBC = ele_stiffness_QBC * area(QBC_i);
% %     mtx.A(mesh.QBC2_elements(1 : 3, i), mesh.QBC2_elements(1 : 3, i)) = ... 
% %         fem.stiffness(mesh.QBC2_elements(1 : 3, i), mesh.QBC2_elements(1 : 3, i)) + ele_stiffness_QBC;
% % 
% %     % add QBC for mass matrix
% %     ele_mass = zeros(3, 3);
% %     for j = 1 : length(qwgts) % quadrature loop
% %         r = rspts(j, 1); % r coordinate
% %         s = rspts(j, 2); % s coordinate
% %         x = mesh.nodes(1, mesh.QBC2_elements(1 : 3, i));
% %         y = mesh.nodes(2, mesh.QBC2_elements(1 : 3, i));
% %         [S, dSdx, dSdy, detJ] = isoparametricMap(x', y', r, s, @P1); % map
% %         wxarea = qwgts(j) * detJ / 2; % weight times det(J)
% %         tensor11_QBC = parameters.math.tensor11_fun(mesh.QBC2_elements(4, i));
% %         ele_mass = ele_mass + (S * S') * (tensor11_QBC * wxarea); % compute and add integrand to MK
% %     end
% %     mtx.B(mesh.QBC2_elements(1 : 3, i), mesh.QBC2_elements(1 : 3, i)) = ... 
% %         mtx.B(mesh.QBC2_elements(1 : 3, i), mesh.QBC2_elements(1 : 3, i)) + ele_mass;
% % end

end