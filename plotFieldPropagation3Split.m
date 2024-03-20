function [P_L, P_R] = plotFieldPropagation3Split(result_CI, mtx, mesh, parameters)

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
result_CI.ratio_qL = zeros(p - 1, 1);
result_CI.ratio_qR = zeros(p - 1, 1);
for i = p - 1 : -1 : 1
    Q_L(:, i) = P_L \ (NLmu' * Q_L(:, i + 1));
    Q_R(:, p - i + 1) = P_R \ (NRmu * Q_R(:, p - i));
    result_CI.ratio_qL(p - i) = norm(Q_L(:, i), 1) / norm(Q_L(:, i + 1), 1);
    result_CI.ratio_qR(p - i) = norm(Q_R(:, p - i + 1), 1) / norm(Q_R(:, p - i), 1);
end
fprintf('mean ratio_qL: %.3f, max eigenvalue: %.3f.\n', mean(result_CI.ratio_qL), max(abs(eigs(P_L \ NLmu'))));
fprintf('mean ratio_qR: %.3f, max eigenvalue: %.3f.\n', mean(result_CI.ratio_qR), max(abs(eigs(P_R \ NRmu))));



if 1

    % split into five parts and rotate
    pp = 8;
    figure
    nodes_origin1 = nodes_origin1 - p * parameters.lattice.shift1;
    for i = 1 : pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, i);
        q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
        abs_q_L = abs(q_L);
        nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
        patch('Faces', elements1(1 : 3, :).', 'Vertices', nodes1.', 'FaceVertexCData', abs_q_L, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin1 = nodes_origin1 + 1.5 * parameters.lattice.a1;
    for i = 1 : pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, pp + i);
        q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
        abs_q_L = abs(q_L);
        nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
        patch('Faces', elements1(1 : 3, :).', 'Vertices', nodes1.', 'FaceVertexCData', abs_q_L, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin1 = nodes_origin1 + 1.5 * parameters.lattice.a1;
    nodes_origin1 = nodes_origin1 - parameters.lattice.shift1;
    for i = 1 : p - 2 * pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, 2 * pp + i);
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
    nodes2 = nodes_origin2 - 2 * pp * parameters.lattice.shift1 + 1.5 * 2 * parameters.lattice.a1;
    nodes2 = nodes2 - parameters.lattice.shift1;
    patch('Faces', elements2(1 : 3, :).', 'Vertices', nodes2.', 'FaceVertexCData', abs_q_C, ...
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on

    nodes_origin3 = nodes_origin3 - p * parameters.lattice.shift1 + 1.5 * 2 * parameters.lattice.a1;
    for i = 1 : p - 2 * pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (p - 2 * pp + 1.5 - 1) * parameters.lattice.shift1;
        nodes3 = nodes3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', nodes3.', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin3 = nodes_origin3 + 1.5 * parameters.lattice.a1;
    for i = 1 : pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, p - 2 * pp + i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', nodes3.', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin3 = nodes_origin3 + 1.5 * parameters.lattice.a1;
    for i = 1 : pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, p - pp + i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', nodes3.', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    colorbar
    axis equal
    axis off

else

    % split into five parts and rotate
    G = givens(parameters.lattice.a2(1), parameters.lattice.a2(2));

    pp = 8;
    blank = 1.8;
    len_arrow = 0.5 * [1, 0];
    start = pp * parameters.lattice.a1;
    text(start(1), start(2), 'arrow');


    figure
    nodes_origin1 = nodes_origin1 - p * parameters.lattice.shift1;
    for i = 1 : pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, i);
        q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
        abs_q_L = abs(q_L);
        nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
        patch('Faces', elements1(1 : 3, :).', 'Vertices', (G * nodes1).', 'FaceVertexCData', abs_q_L, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin1 = nodes_origin1 + blank * parameters.lattice.a1 - 0 * parameters.lattice.shift1;
    for i = 1 : pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, pp + i);
        q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
        abs_q_L = abs(q_L);
        nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
        patch('Faces', elements1(1 : 3, :).', 'Vertices', (G * nodes1).', 'FaceVertexCData', abs_q_L, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin1 = nodes_origin1 + blank * parameters.lattice.a1;
    nodes_origin1 = nodes_origin1 - parameters.lattice.shift1 + 1 * parameters.lattice.shift1;
    for i = 1 : p - 2 * pp
        q_L = zeros(n1, 1);
        q_L(index1.free) = Q_L(:, 2 * pp + i);
        q_L(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_L(index1.node_quasiperiodic1);
        abs_q_L = abs(q_L);
        nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
        patch('Faces', elements1(1 : 3, :).', 'Vertices', (G * nodes1).', 'FaceVertexCData', abs_q_L, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    q_C = zeros(n2, 1);
    q_C(index2.free) = qC_inf;
    q_C(index2.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_C(index2.node_quasiperiodic1);
    abs_q_C = abs(q_C);
    nodes2 = nodes_origin2 - 2 * pp * parameters.lattice.shift1 + blank * 2 * parameters.lattice.a1;
    nodes2 = nodes2 - parameters.lattice.shift1 + 1 * parameters.lattice.shift1;
    patch('Faces', elements2(1 : 3, :).', 'Vertices', (G * nodes2).', 'FaceVertexCData', abs_q_C, ...
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on

    nodes_origin3 = nodes_origin3 - p * parameters.lattice.shift1 + blank * 2 * parameters.lattice.a1 + 0 * parameters.lattice.shift1;
    for i = 1 : p - 2 * pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (p - 2 * pp + 2.5 - 1) * parameters.lattice.shift1;
        nodes3 = nodes3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', (G * nodes3).', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin3 = nodes_origin3 + blank * parameters.lattice.a1 - 0 * parameters.lattice.shift1;
    for i = 1 : pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, p - 2 * pp + i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', (G * nodes3).', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin3 = nodes_origin3 + blank * parameters.lattice.a1 - 0 * parameters.lattice.shift1;
    for i = 1 : pp
        q_R = zeros(n3, 1);
        q_R(index3.free) = Q_R(:, p - pp + i);
        q_R(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index3.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
        patch('Faces', elements3(1 : 3, :).', 'Vertices', (G * nodes3).', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    colorbar
    axis equal
    axis off

end

end