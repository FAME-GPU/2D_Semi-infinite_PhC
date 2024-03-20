function [P_L, P_R] = plotFieldPropagation3(result_CI, mtx, mesh, parameters)

gamma = parameters.gamma;
tol = parameters.SDA.tol;
K_C = mtx.K_C;
M_C = mtx.M_C;
K_L = mtx.K_L;
M_L = mtx.M_L;
E_L = mtx.E_L;
F_L = mtx.F_L;
E_LB = mtx.E_LB;
F_LB = mtx.F_LB;
K_R = mtx.K_R;
M_R = mtx.M_R;
E_R = mtx.E_R;
F_R = mtx.F_R;
E_RB = mtx.E_RB;
F_RB = mtx.F_RB;

%% plot non-referenced fields
% plot q1_inf
mu = result_CI.mu;
qC_inf = result_CI.qC_inf;
p = parameters.num_propagation;
nodes_origin1 = mesh.geometry1.nodes;
index1 = mesh.geometry1.index;
elements1 = mesh.geometry1.elements;
n1 = mesh.geometry1.n_nodes;
m1 = mesh.geometry1.n_free_nodes;
nodes_origin2 = mesh.geometry2.nodes;
index2 = mesh.geometry2.index;
elements2 = mesh.geometry2.elements;
n2 = mesh.geometry2.n_nodes;
m2 = mesh.geometry2.n_free_nodes;
nodes_origin3 = mesh.geometry3.nodes;
index3 = mesh.geometry3.index;
elements3 = mesh.geometry3.elements;
n3 = mesh.geometry3.n_nodes;
m3 = mesh.geometry3.n_free_nodes;

% plot q1_inf propagation
% ACmu = mu * M_C - K_C;
ALmu = mu * M_L - K_L;
ARmu = mu * M_R - K_R;
NLmu = E_L - mu * F_L;
NRmu = E_R - mu * F_R;
NLBmu = E_LB - mu * F_LB;
NRBmu = E_RB - mu * F_RB;
% [invGLmu, ~] = FPI(NLmu', ALmu, 0.25, tol);
% [invGRmu, ~] = FPI(NRmu, ARmu, 0.25, tol);
[invGLmu, ~] = SDA_CI(NLmu', ALmu, tol);
[invGRmu, ~] = SDA_CI(NRmu, ARmu, tol);
P_L = ALmu - NLmu * (invGLmu \ NLmu');
P_R = ARmu - NRmu' * (invGRmu \ NRmu);
Q_L = zeros(m1, p);
Q_R = zeros(m3, p);
Q_L(:, p) = P_L \ (NLBmu' * qC_inf);
Q_R(:, 1) = P_R \ (NRBmu * qC_inf);
for i = p - 1 : -1 : 1
    Q_L(:, i) = P_L \ (NLmu' * Q_L(:, i + 1));
    Q_R(:, p - i + 1) = P_R \ (NRmu * Q_R(:, p - i));
end

figure
nodes_origin1 = nodes_origin1 - p * parameters.lattice.shift1;
for i = 1 : p
    q_L = zeros(n1, 1);
    q_L(index1.free) = Q_L(:, i);
    q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
    abs_q_L = abs(q_L);
    nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
    patch('Faces', elements1(1 : 3, :).', 'Vertices', nodes1.', 'FaceVertexCData', abs_q_L, ... 
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end

q_C = zeros(n2, 1);
q_C(index2.free) = qC_inf;
q_C(index2.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_C(index2.node_quasiperiodic1);
abs_q_C = abs(q_C);
nodes2 = nodes_origin2; % parameters.lattice.shift;
patch('Faces', elements2(1 : 3, :).', 'Vertices', nodes2.', 'FaceVertexCData', abs_q_C, ... 
    'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on

nodes_origin3 = nodes_origin3 + 1.5 * parameters.lattice.shift1;
for i = 1 : p
    q_R = zeros(n3, 1);
    q_R(index3.free) = Q_R(:, i);
    q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
    abs_q_R = abs(q_R);
    nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
    patch('Faces', elements3(1 : 3, :).', 'Vertices', nodes3.', 'FaceVertexCData', abs_q_R, ... 
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end
colorbar
axis equal
axis off


% %% plot referenced fields
% % plot q1_inf
% figure
% qC_inf = result_CI.q1_inf_reference;
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
% [invG11_mu, ~] = FPI(EmF, mu * mtx.M_R - mtx.K_R, 0.25, parameters.SDA.tol);
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
% if 1
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


