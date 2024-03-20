function [result_CI, result_select] = computeSelectedEigenpairs(mtx, result, parameters)
% Use contour integral to compute the eigenpairs of infinite structures

K_C = mtx.K_C;
M_C = mtx.M_C;
K_R = mtx.K_R;
M_R = mtx.M_R;
E_R = mtx.E_R;
F_R = mtx.F_R;
E_RB = mtx.E_RB;
F_RB = mtx.F_RB;

m = size(K_C, 1);
num_points = parameters.contour_integral.num_points;
radius = parameters.contour_integral.radius;
tol = parameters.SDA.tol;

ew_finite_select = result{parameters.selected_index_wave_vec}.ew_finite(parameters.selected_index_eigencurve);
ev_finite_select = result{parameters.selected_index_wave_vec}.ev_finite(:, parameters.selected_index_eigencurve);

%% contour integral
time_contour_integral = tic;
ev = 0;
ew_ev = 0;
time_SDA = zeros(num_points, 1);
time_cSDA = zeros(num_points, 1);
step_SDA = zeros(num_points, 1);
step_cSDA = zeros(num_points, 1);
for i = 1 : num_points
    z = ew_finite_select + radius * exp(1i * 2 * pi * i / num_points);
    ACz = z * M_C - K_C;
    ARz = z * M_R - K_R;
    NRz = E_R - z * F_R;
    NRBz = E_RB - z * F_RB;

    tmp_time_SDA_R = tic;
    [invGRz, ~, step_SDA(i)] = SDA_CI(NRz, ARz, tol);
    invGCz = ACz - NRBz' * (invGRz \ NRBz);
    time_SDA(i) = toc(tmp_time_SDA_R);

    % % test SDA_cyclic
    % p = size(K_R, 1) / 2;
    % permute = symrcm(ARz);
    % new_ARz = ARz(permute, permute);
    % new_ARz2 = new_ARz(1 : p, 1 : p);
    % new_ARz1 = new_ARz(p + 1 : end, p + 1 : end);
    % new_NRz1 = new_ARz(p + 1 : end, 1 : p);
    % new_NRz2 = NRz(1 : p, p + 1 : end);
    % tmp_time_cSDA = tic;
    % A = new_NRz2 * (new_ARz2 \ new_NRz1);
    % Bs = new_NRz1' * (new_ARz2 \ new_NRz2);
    % Q = new_ARz1 - new_NRz1' * (new_ARz2 \ new_NRz1);
    % P = new_NRz2 * (new_ARz2 \ new_NRz2');
    % [~, ~, invGRz1, ~, step_cSDA(i)] = SDA_cyclic_CI(A, Bs, Q, P, 1e-10);
    % time_cSDA(i) = toc(tmp_time_cSDA);
    % Lz1 = new_NRz2' * (invGRz1 \ new_NRz2);
    % % invGRz2 = new_ARz2 - Lz1;
    % % Lz2 = new_NRz1' * (invGRz2 \ new_NRz1);
    % 
    % new_ACz = ACz(permute, permute);
    % new_ACz0 = new_ACz(1 : p, 1 : p);
    % invGCz0 = new_ACz0 - Lz1;
    
    

    dt = 2 * pi / num_points;
    dz = 1i * radius * exp(1i * 2 * pi * i / num_points) * dt;
    ele_integral = invGCz \ (dz * speye(m));
    ev = ev + ele_integral;
    ew_ev = ew_ev + z * ele_integral;
end
time_SDA = sort(time_SDA);
time_SDA = time_SDA(1 : end - 3);
time_cSDA = sort(time_cSDA);
time_cSDA = time_cSDA(1 : end - 3);
fprintf('%d %.1e\n', mode(step_SDA), mean(time_SDA));
fprintf('%d %.1e %.1f\n', mode(step_cSDA), mean(time_cSDA), mean(time_SDA) / mean(time_cSDA));
ev = ev / (2 * pi * 1i);
ew_ev = ew_ev / (2 * pi * 1i);
time_contour_integral = toc(time_contour_integral);

ew_one = eig(ev);
[~, index] = sort(abs(ew_one), 'descend'); 
ew_one = ew_one(index);
isRank1 = 1 - abs(ew_one(2) / ew_one(1));
fprintf('Is rank 1: %.12f.\n', isRank1);

ew_mu = eig(ew_ev);
[~, index] = sort(abs(ew_mu), 'descend'); 
ew_mu = ew_mu(index);
isRank1_2 = 1 - abs(ew_mu(2) / ew_mu(1));
fprintf('Is rank 1_2: %.12f.\n', isRank1_2);

mu = real(ew_mu(1) / ew_one(1));

[qC, ~] = qr(ev);
qC = qC(:, 1);

qC_finite = ev_finite_select(1 : m);
qC_finite = qC_finite / norm(qC_finite);
ratio_qC = qC_finite ./ qC;
isEigenvectorC = 1 - std(ratio_qC);
fprintf('Is eigenvector: %.12f.\n', isEigenvectorC);

result_CI.ev_inf = ev;
result_CI.ew_ev_inf = ew_ev;
result_CI.qC_inf = qC;
result_CI.ratio_qC = ratio_qC;
result_CI.isRank1 = isRank1;
result_CI.isEigenvectorC = isEigenvectorC;
result_CI.time_contour_integral = time_contour_integral;
result_CI.time_SDA = time_SDA;
result_CI.ew_one = ew_one;
result_CI.ew_mu = ew_mu;
result_CI.mu = mu;
result_CI.mu_normalized = sqrt(mu) * parameters.lattice.a / (2 * pi);

result_select.qC_finite = qC_finite;
result_select.ev_finite_select = ev_finite_select;
result_select.ew_finite_select = ew_finite_select;

end

% function [result_CI, result_select] = computeSelectedEigenpairs(mtx, result, parameters)
% % Use contour integral to compute the eigenpairs of infinite structures
% 
% K_C = mtx.K_C;
% M_C = mtx.M_C;
% K_R = mtx.K_R;
% E_R = mtx.E_R;
% M_R = mtx.M_R;
% F_R = mtx.F_R;
% 
% n = size(K_R, 1);
% num_points = parameters.contour_integral.num_points;
% radius = parameters.contour_integral.radius;
% tol = parameters.SDA.tol;
% 
% ew_finite_select = result{parameters.selected_index_wave_vec}.ew_finite(parameters.selected_index_eigencurve);
% ev_finite_select = result{parameters.selected_index_wave_vec}.ev_finite(:, parameters.selected_index_eigencurve);
% 
% %% contour integral
% time_contour_integral = tic;
% ev_inf = 0;
% ew_ev_inf = 0;
% time_SDA = zeros(num_points, 1);
% for i = 1 : num_points
%     z = ew_finite_select + radius * exp(1i * 2 * pi * i / num_points);
%     tmp_time_SDA = tic;
%     [invG11z, ~, step] = SDA_CI(E_R - z * F_R, z * M_R - K_R, tol);
%     time_SDA(i) = toc(tmp_time_SDA);
%     fprintf('Residual of SDA_CI: %.3e with %d steps. Time_SDA: %.3e s.\n', ...
%         norm(invG11z + (E_R - z * F_R)' * (invG11z \ (E_R - z * F_R)) - (z * M_R - K_R), 1), step, time_SDA(i));
%     % G11z = inv(invG11z);
%     invG11z_0 = z * M_C - K_C - (E_R - z * F_R)' * (invG11z \ (E_R - z * F_R));
%     % G11z_0 = inv(invG11z_0);
%     dt = 2 * pi / num_points;
%     dz = 1i * radius * exp(1i * 2 * pi * i / num_points) * dt;
%     ele_integral = invG11z_0 \ (dz * speye(n));
%     ev_inf = ev_inf + ele_integral;
%     ew_ev_inf = ew_ev_inf + z * ele_integral;
% end
% ev_inf = ev_inf / (2 * pi * 1i);
% ew_ev_inf = ew_ev_inf / (2 * pi * 1i);
% time_contour_integral = toc(time_contour_integral);
% 
% ew_one = eig(ev_inf);
% [~, index] = sort(abs(ew_one), 'descend'); 
% ew_one = ew_one(index);
% isRank1 = 1 - abs(ew_one(2) / ew_one(1));
% fprintf('Is rank 1: %.12f.\n', isRank1);
% 
% ew_mu = eig(ew_ev_inf);
% [~, index] = sort(abs(ew_mu), 'descend'); 
% ew_mu = ew_mu(index);
% isRank1_2 = 1 - abs(ew_mu(2) / ew_mu(1));
% fprintf('Is rank 1_2: %.12f.\n', isRank1_2);
% 
% mu = real(ew_mu(1) / ew_one(1));
% fprintf('mu - 0.436 = %.12f.\n', norm(sqrt(mu) / (2 * pi) - 0.436));
% 
% [q1_inf, ~] = qr(ev_inf);
% q1_inf = q1_inf(:, 1);
% 
% q1_finite = ev_finite_select(1 : n);
% q1_finite = q1_finite / norm(q1_finite);
% ratio_q1 = q1_finite ./ q1_inf;
% isEigenvector1 = 1 - std(ratio_q1);
% fprintf('Is eigenvector: %.12f.\n', isEigenvector1);
% 
% result_CI.ev_inf = ev_inf;
% result_CI.ew_ev_inf = ew_ev_inf;
% result_CI.q1_inf = q1_inf;
% result_CI.ratio_q1 = ratio_q1;
% result_CI.isRank1 = isRank1;
% result_CI.isEigenvector1 = isEigenvector1;
% result_CI.time_contour_integral = time_contour_integral;
% result_CI.time_SDA = time_SDA;
% result_CI.mu = mu;
% result_CI.ew_one = ew_one;
% result_CI.ew_mu = ew_mu;
% 
% result_select.q1_finite = q1_finite;
% result_select.ev_finite_select = ev_finite_select;
% result_select.ew_finite_select = ew_finite_select;
% 
% %% construct reference matrices
% num_blocks = parameters.num_blocks;
% 
% A_finite_reference = kron(speye(num_blocks - 1), K_R);
% row_index_E = 2 : num_blocks;
% col_index_E = 1 : num_blocks - 1;
% E_finite_reference = sparse(row_index_E, col_index_E, 1, num_blocks, num_blocks, num_blocks - 1);
% E_finite_reference = kron(E_finite_reference, E_R);
% A_finite_reference = blkdiag(K_C, A_finite_reference);
% A_finite_reference = A_finite_reference + E_finite_reference + E_finite_reference';
% 
% B_finite_reference = kron(speye(num_blocks - 1), M_R);
% row_index_F = 2 : num_blocks;
% col_index_F = 1 : num_blocks - 1;
% F_finite_reference = sparse(row_index_F, col_index_F, 1, num_blocks, num_blocks, num_blocks - 1);
% F_finite_reference = kron(F_finite_reference, F_R);
% B_finite_reference = blkdiag(M_C, B_finite_reference);
% B_finite_reference = B_finite_reference + F_finite_reference + F_finite_reference';
% 
% [ev_finite_reference, ew_finite_reference] = eigs(A_finite_reference, B_finite_reference, ... 
%     parameters.n_eigenvalues, 'smallestabs', 'Tolerance', 1e-12);
% [ew_finite_reference, index] = sort(real(diag(ew_finite_reference)), "ascend");
% ev_finite_reference = ev_finite_reference(:, index);
% ev_finite_reference_select = ev_finite_reference(:, parameters.selected_index_eigencurve);
% ew_finite_reference_select = ew_finite_reference(parameters.selected_index_eigencurve);
% fprintf('|mu - ew_finite_reference_select| = %.12f.\n', abs(mu - ew_finite_reference_select));
% 
% %% contour integral for reference eigenvalue
% time_contour_integral_reference = tic;
% ev_inf_reference = 0;
% ew_ev_inf_reference = 0;
% for i = 1 : num_points
%     z = ew_finite_reference_select + radius * exp(1i * 2 * pi * i / num_points);
%     time_SDA = tic;
%     [invG11z, ~, step] = SDA_CI(E_R - z * F_R, z * M_R - K_R, tol);
%     time_SDA = toc(time_SDA);
%     fprintf('Residual of SDA_CI: %.3e with %d steps. Time_SDA: %.3e s.\n', ...
%         norm(invG11z + (E_R - z * F_R)' * (invG11z \ (E_R - z * F_R)) - (z * M_R - K_R), 1), step, time_SDA);
%     % G11z = inv(invG11z);
%     invG11z_0 = z * M_C - K_C - (E_R - z * F_R)' * (invG11z \ (E_R - z * F_R));
%     % G11z_0 = inv(invG11z_0);
%     dt = 2 * pi / num_points;
%     dz = 1i * radius * exp(1i * 2 * pi * i / num_points) * dt;
%     ele_integral = invG11z_0 \ (dz * speye(n));
%     ev_inf_reference = ev_inf_reference + ele_integral;
%     ew_ev_inf_reference = ew_ev_inf_reference + z * ele_integral;
% end
% ev_inf_reference = ev_inf_reference / (2 * pi * 1i);
% ew_ev_inf_reference = ew_ev_inf_reference / (2 * pi * 1i);
% time_contour_integral_reference = toc(time_contour_integral_reference);
% 
% ew_one_reference = eig(ev_inf_reference);
% [~, index] = sort(abs(ew_one_reference), 'descend'); 
% ew_one_reference = ew_one_reference(index);
% isRank1_reference = 1 - abs(ew_one_reference(2) / ew_one_reference(1));
% fprintf('Is rank 1 reference: %.12f.\n', isRank1_reference);
% 
% ew_mu_reference = eig(ew_ev_inf_reference);
% [~, index] = sort(abs(ew_mu_reference), 'descend'); 
% ew_mu_reference = ew_mu_reference(index);
% isRank1_2_reference = 1 - abs(ew_mu_reference(2) / ew_mu_reference(1));
% fprintf('Is rank 1_2 reference: %.12f.\n', isRank1_2_reference);
% 
% mu_reference = real(ew_mu_reference(1) / ew_one_reference(1));
% % fprintf('mu_reference - 0.436 = %.12f.\n', norm(sqrt(mu_reference) / (2 * pi) - 0.436));
% 
% [q1_inf_reference, ~] = qr(ev_inf_reference);
% q1_inf_reference = q1_inf_reference(:, 1);
% 
% q1_finite_reference = ev_finite_reference_select(1 : n);
% q1_finite_reference = q1_finite_reference / norm(q1_finite_reference);
% ratio_q1_reference = q1_finite_reference ./ q1_inf_reference;
% isEigenvector1_reference = 1 - std(ratio_q1_reference);
% fprintf('Is eigenvector: %.12f.\n', isEigenvector1_reference);
% 
% result_CI.ev_inf_reference = ev_inf_reference;
% result_CI.ew_ev_inf_reference = ew_ev_inf_reference;
% result_CI.q1_inf_reference = q1_inf_reference;
% result_CI.ratio_q1_reference = ratio_q1_reference;
% result_CI.isRank1_reference = isRank1_reference;
% result_CI.isEigenvector1_reference = isEigenvector1_reference;
% result_CI.time_contour_integral_reference = time_contour_integral_reference;
% result_CI.mu_reference = mu_reference;
% result_CI.ew_one_reference = ew_one_reference;
% result_CI.ew_mu_reference = ew_mu_reference;
% 
% result_select.q1_finite_reference = q1_finite_reference;
% result_select.ev_finite_reference_select = ev_finite_reference_select;
% result_select.ew_finite_reference_select = ew_finite_reference_select;
% 
% % %% plot eigenvectors
% % free_nodes = mesh.nodes(:, mesh.new.index_free);
% % 
% % figure
% % scatter(free_nodes(1, :), free_nodes(2, :), length(free_nodes(1, :)), abs(q1_finite), 'filled', 'Marker', 'square');
% % colormap('default')
% % colorbar
% % axis equal
% % axis off
% % grid on
% % 
% % figure
% % scatter(free_nodes(1, :), free_nodes(2, :), length(free_nodes(1, :)), abs(q1_inf), 'filled', 'Marker', 'square')
% % colormap('default')
% % colorbar
% % axis equal
% % axis off
% % grid on
% 
% end