function eigenpairs = computeEigenpairs(mtx, parameters)
% function eigenpairs = computeEigenpairs(mtx, mesh, reference, parameters)
% This function computes the eigenpairs of supercell structures


% num_points = parameters.contour_integral.num_points;
% radius = parameters.contour_integral.radius;
% tol = parameters.SDA.tol;
% n = size(mtx.A0, 1);
% num_blocks = parameters.num_blocks;
% ind_select = 23;  

n = size(mtx.KK, 1);
n_eigenvalues = parameters.n_eigenvalues;

% A0 = mtx.A0;
% A = mtx.A;
% Ae = mtx.Ae;
% E0 = mtx.E0;
% E = mtx.E;
% Ee = mtx.Ee;
% B0 = mtx.B0;
% B = mtx.B;
% Be = mtx.Be;
% F0 = mtx.F0;
% F = mtx.F;
% Fe = mtx.Fe;

% A_finite = sparse(n * num_blocks, n * num_blocks);
% for jj = 1 : num_blocks - 1
%     A_finite((jj - 1) * n + 1 : jj * n, (jj - 1) * n + 1 : jj * n) = A;
%     A_finite((jj - 1) * n + 1 : jj * n, jj * n + 1 : (jj + 1) * n) = E';
%     A_finite(jj * n + 1 : (jj + 1) * n, (jj - 1) * n + 1 : jj * n) = E;
% end
% A_finite((num_blocks - 1) * n + 1 : num_blocks * n, (num_blocks - 1) * n + 1 : num_blocks * n) = Ae;
% A_finite = blkdiag(A0, A_finite);
% A_finite(1 : n, n + 1 : 2 * n) = E';
% A_finite(n + 1 : 2 * n, 1 : n) = E;
% 
% B_finite = sparse(n * num_blocks, n * num_blocks);
% for jj = 1 : num_blocks - 1
%     B_finite((jj - 1) * n + 1 : jj * n, (jj - 1) * n + 1 : jj * n) = B;
%     B_finite((jj - 1) * n + 1 : jj * n, jj * n + 1 : (jj + 1) * n) = F';
%     B_finite(jj * n + 1 : (jj + 1) * n, (jj - 1) * n + 1 : jj * n) = F;
% end
% B_finite((num_blocks - 1) * n + 1 : num_blocks * n, (num_blocks - 1) * n + 1 : num_blocks * n) = Be;
% B_finite = blkdiag(B0, B_finite);
% B_finite(1 : n, n + 1 : 2 * n) = F';
% B_finite(n + 1 : 2 * n, 1 : n) = F;

% A_finite = kron(speye(num_blocks - 2), A);
% A_finite = blkdiag(A0, A_finite, Ae);
% row_index_E = 2 : num_blocks - 2;
% col_index_E = 1 : num_blocks - 3;
% E_finite = sparse(row_index_E, col_index_E, 1, num_blocks - 2, num_blocks - 2, num_blocks - 3);
% E_finite = kron(E_finite, E);
% E_finite = blkdiag(sparse(size(A0, 1), size(A0, 2)), E_finite, sparse(size(Ae, 1), size(Ae, 2)));
% E_finite(size(A0, 1) + 1 : size(A0, 1) + size(A, 1), 1 : size(A0, 2)) = E0;
% E_finite(end - size(Ae, 1) + 1 : end, end - size(Ae, 2) - size(A, 2) + 1 : end - size(Ae, 2)) = Ee;
% A_finite = A_finite + E_finite + E_finite';
% 
% B_finite = kron(speye(num_blocks - 2), B);
% B_finite = blkdiag(B0, B_finite, Be);
% row_index_F = 2 : num_blocks - 2;
% col_index_F = 1 : num_blocks - 3;
% F_finite = sparse(row_index_F, col_index_F, 1, num_blocks - 2, num_blocks - 2, num_blocks - 3);
% F_finite = kron(F_finite, F);
% F_finite = blkdiag(sparse(size(B0, 1), size(B0, 2)), F_finite, sparse(size(Be, 1), size(Be, 2)));
% F_finite(size(B0, 1) + 1 : size(B0, 1) + size(B, 1), 1 : size(B0, 2)) = F0;
% F_finite(end - size(Be, 1) + 1 : end, end - size(Be, 2) - size(B, 2) + 1 : end - size(Be, 2)) = Fe;
% B_finite = B_finite + F_finite + F_finite';


time_eigs = tic;
[ev_finite, ew_finite] = eigs(mtx.KK, mtx.MM, n_eigenvalues, 'smallestabs', 'Tolerance', 1e-12);
time_eigs = toc(time_eigs);
[ew_finite, index] = sort(diag(real(ew_finite)), 'ascend');
ev_finite = ev_finite(:, index);

% residual = zeros(n_eigenvalues, 1);
% for i = 1 : n_eigenvalues
%     residual(i) = norm(A_finite * ev_finite(:, i) - ew_finite(i) * (B_finite * ev_finite(:, i)), inf);
% end

eigenpairs.ew_finite = ew_finite;
eigenpairs.ev_finite = ev_finite;
eigenpairs.time_eigs = time_eigs;

% figure
% plot(1 : n_eigenvalues, real(ew_finite), 'bo')
% message = "Input the selected index: ";
% ind_select = input(message);
% ew_select = ew_finite(ind_select);

% %% contour integral
% tic
% ev_inf = 0;
% % ew_ev_CI = 0;
% for i = 1 : num_points
%     z = ew_select + radius * exp(1i * 2 * pi * i / num_points);
%     time_SDA = tic;
%     [invG11z, ~, step] = SDA_CI(E - z * F, z * B - A, tol);
%     time_SDA = toc(time_SDA);
%     fprintf('Residual of SDA_CI: %.3e with %d steps. Time_SDA: %.3e s.\n', ...
%         norm(invG11z + (E - z * F)' * (invG11z \ (E - z * F)) - (z * B - A), 1), step, time_SDA);
%     % G11z = inv(invG11z);
%     invG11z_0 = z * B0 - A0 - (E - z * F)' * (invG11z \ (E - z * F));
%     % G11z_0 = inv(invG11z_0);
%     dt = 2 * pi / num_points;
%     dz = 1i * radius * exp(1i * 2 * pi * i / num_points) * dt;
%     ev_inf = ev_inf + invG11z_0 \ (dz * speye(n));
%     % ew_ev_CI = ew_ev_CI + z * invG11z_0 \ (dz * speye(n));
% end
% ev_inf = ev_inf / (2 * pi * 1i);
% % ew_ev_CI = ew_ev_CI / (2 * pi * 1i);
% toc
% 
% ew_rank1 = eig(ev_inf);
% [~, index] = sort(abs(ew_rank1), 'descend'); 
% ew_rank1 = ew_rank1(index);
% isRank1 = 1 - abs(ew_rank1(2) / ew_rank1(1));
% fprintf('Is rank 1: %.12f.\n', isRank1);
% 
% [q1_inf, ~] = qr(ev_inf);
% q1_inf = q1_inf(:, 1);
% 
% ev_finite_select = ev_finite(:, ind_select);
% q1_finite = ev_finite_select(1 : n);
% ratio = q1_finite ./ q1_inf;
% isEigenvector = 1 - (norm(ratio, 1) / n - norm(ratio, inf)) / norm(ratio, inf);
% fprintf('Is eigenvector: %.12f.\n', isEigenvector);


% %% plot eigenvectors
% free_nodes = mesh.nodes(:, mesh.new.index_free);
% figure
% scatter(free_nodes(1, :), free_nodes(2, :), length(free_nodes(1, :)), abs(true_q1)', 'filled')
% colormap('default')
% colorbar
% axis equal
% axis off
% grid on
% 
% figure
% scatter(free_nodes(1, :), free_nodes(2, :), length(free_nodes(1, :)), abs(q1)', 'filled')
% colormap('default')
% colorbar
% axis equal
% axis off
% grid on
% 
% figure
% for i = 1 : num_blocks
%     q1_i = abs(ev_select((i - 1) * n + 1 : i * n));
%     nodes_x = free_nodes(1, :) + (i - 1) * parameters.lattice.shift(1);
%     nodes_y = free_nodes(2, :) + (i - 1) * parameters.lattice.shift(2);
%     scatter(nodes_x, nodes_y, n, q1_i, 'filled'), hold on
% end
% colorbar
% axis equal
% axis off
% grid on


end


