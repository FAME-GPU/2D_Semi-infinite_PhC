function [propagate_matrix, EmF] = plotFieldPropagation(result_CI, mtx, mesh, parameters)

%% plot non-referenced fields
% plot q1_inf
figure
qC_inf = result_CI.qC_inf;
nodes_origin = mesh.geometry1.nodes;
index = mesh.geometry1.index;
elements = mesh.geometry1.elements;
n = mesh.geometry1.n_nodes;
m = mesh.geometry1.n_free_nodes;

% nodes_origin = mesh.nodes;
% n = mesh.n_nodes;
% m = mesh.n_free_nodes;
q1_inf_all_nodes = zeros(n, 1);
q1_inf_all_nodes(index.free) = qC_inf;
q1_inf_all_nodes(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_inf_all_nodes(index.node_quasiperiodic1);
abs_q1_inf_all_nodes = abs(q1_inf_all_nodes);
patch('Faces', elements(1 : 3, :).', 'Vertices', nodes_origin.', 'FaceVertexCData', abs_q1_inf_all_nodes, ...
    'FaceColor', 'interp', 'EdgeColor', 'interp');
colorbar
axis equal
axis off

% plot q1_inf propagation
p = parameters.num_propagation;
mu = result_CI.mu;
EmF = mtx.E_R - mu * mtx.F_R;
% [invG11_mu, ~] = FPI(EmF, mu * mtx.M_R - mtx.K_R, 0.25, parameters.SDA.tol);
[invG11_mu, ~] = SDA_CI(EmF, mu * mtx.M_R - mtx.K_R, parameters.SDA.tol);
propagate_matrix = mu * mtx.M_R - mtx.K_R - EmF' * (invG11_mu \ EmF);
Q_inf = zeros(m, p);
Q_inf(:, 1) = qC_inf;
for i = 1 : p - 1
    Q_inf(:, i + 1) = propagate_matrix \ (EmF * Q_inf(:, i));
end
% plot
figure
for i = 1 : p
    q_inf_all_nodes = zeros(n, 1);
    q_inf_all_nodes(index.free) = Q_inf(:, i);
    q_inf_all_nodes(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_inf_all_nodes(index.node_quasiperiodic1);
    abs_q_inf_all_nodes = abs(q_inf_all_nodes);
    nodes = nodes_origin + (i - 1) * parameters.lattice.shift;
    patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q_inf_all_nodes, ...
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end
colorbar
axis equal
axis off


% %% plot referenced fields
% % plot q1_inf
% figure
% qC_inf = result_CI.qC_inf_reference;
% q1_inf_all_nodes = zeros(n, 1);
% q1_inf_all_nodes(index.free) = qC_inf;
% q1_inf_all_nodes(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_inf_all_nodes(index.node_quasiperiodic1);
% abs_q1_inf_all_nodes = abs(q1_inf_all_nodes);
% patch('Faces', elements(1 : 3, :).', 'Vertices', nodes_origin.', 'FaceVertexCData', abs_q1_inf_all_nodes, ...
%     'FaceColor', 'interp', 'EdgeColor', 'interp');
% colorbar
% axis equal
% axis off
%
% % plot q1_inf propagation
% p = parameters.num_propagation;
% mu = result_CI.mu_reference;
% EmF = mtx.E_R - mu * mtx.F_R;
% % [invG11_mu, ~] = FPI(EmF, mu * mtx.M_R - mtx.K_R, 0.25, parameters.SDA.tol);
% [invG11_mu, ~] = SDA_CI(EmF, mu * mtx.M_R - mtx.K_R, parameters.SDA.tol);
% propagate_matrix = mu * mtx.M_R - mtx.K_R - EmF' * (invG11_mu \ EmF);
% Q_inf = zeros(m, p);
% Q_inf(:, 1) = qC_inf;
% for i = 1 : p - 1
%     Q_inf(:, i + 1) = propagate_matrix \ (EmF * Q_inf(:, i));
% end
% % plot
% figure
% for i = 1 : p
%     q_inf_all_nodes = zeros(n, 1);
%     q_inf_all_nodes(index.free) = Q_inf(:, i);
%     q_inf_all_nodes(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_inf_all_nodes(index.node_quasiperiodic1);
%     abs_q_inf_all_nodes = abs(q_inf_all_nodes);
%     nodes = nodes_origin + (i - 1) * parameters.lattice.shift;
%     patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q_inf_all_nodes, ...
%         'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
% end
% colorbar
% axis equal
% axis off
% 
% if 0
%     K = propagate_matrix \ EmF;
%     figure
%     [Psi, ew_K] = eig(full(K));
%     [~, ind] = sort(abs(diag(ew_K)), 'descend');
%     Psi = Psi(:, ind);
%     q1_all = zeros(n, 1);
%     q1 = Psi(1 : m, 1);
%     q1 = q1 / norm(q1);
%     q1_all(index.free) = q1;
%     q1_all(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_all(index.node_quasiperiodic1);
%     abs_q1_all = abs(q1_all);
%     nodes = nodes_origin;
%     patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q1_all, ...
%         'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
%     colorbar
%     axis equal
%     axis off
% 
%     % plot y_i_inf
%     tmp = zeros(parameters.num_propagation, 1);
%     for i = 1 : parameters.num_propagation
%         tmp(i) = (Q_inf(:, i)' * q1) / (norm(Q_inf(:, i)) * norm(q1));
%     end
% 
%     v = VideoWriter('DecayOfYk', 'MPEG-4');
%     open(v);
%     figure
%     for i = 1 : parameters.num_propagation
%         clf;
%         qC_inf = Q_inf(:, i);
%         q1_inf_all_nodes = zeros(n, 1);
%         q1_inf_all_nodes(index.free) = qC_inf;
%         q1_inf_all_nodes(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_inf_all_nodes(index.node_quasiperiodic1);
%         abs_q1_inf_all_nodes = abs(q1_inf_all_nodes);
%         patch('Faces', elements(1 : 3, :).', 'Vertices', nodes_origin.', 'FaceVertexCData', abs_q1_inf_all_nodes, ...
%             'FaceColor', 'interp', 'EdgeColor', 'interp');
%         colorbar
%         axis equal
%         axis off
% 
%         tmp = getframe(gcf);
%         frame = tmp.cdata;
%         writeVideo(v, frame);
%         pause(0.1);
%     end
% end

end


