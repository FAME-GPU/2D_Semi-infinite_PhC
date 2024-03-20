function P_R = plotFieldPropagationSplit(result_CI, mtx, mesh, parameters)

gamma = parameters.gamma;
tol = parameters.SDA.tol;
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
nodes_origin = mesh.geometry1.nodes;
index = mesh.geometry1.index;
elements = mesh.geometry1.elements;
n = mesh.geometry1.n_nodes;
m = mesh.geometry1.n_free_nodes;

% plot q1_inf propagation
% ACmu = mu * M_C - K_C;
ARmu = mu * M_R - K_R;
NRmu = E_R - mu * F_R;
[invGRmu, ~] = SDA_CI(NRmu, ARmu, tol);
P_R = ARmu - NRmu' * (invGRmu \ NRmu);
Q_R = zeros(m, p);
Q_R(:, 1) = qC_inf;
result_CI.ratio_q = zeros(p - 1, 1);
for i = p - 1 : -1 : 1
    Q_R(:, p - i + 1) = P_R \ (NRmu * Q_R(:, p - i));
    result_CI.ratio_q(p - i) = norm(Q_R(:, p - i + 1), 1) / norm(Q_R(:, p - i), 1);
end
fprintf('mean ratio_q: %.3f, max eigenvalue: %.3f.\n', mean(result_CI.ratio_q), max(abs(eigs(P_R \ NRmu))));


if 1

    % split into five parts and rotate
    pp = 10;
    figure
    for i = 1 : pp
        q_R = zeros(n, 1);
        q_R(index.free) = Q_R(:, i);
        q_R(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes = nodes_origin + (i - 1) * parameters.lattice.shift; % parameters.lattice.shift;
        patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    nodes_origin = nodes_origin + 2.5 * parameters.lattice.a1;
    for i = 1 : pp
        q_R = zeros(n, 1);
        q_R(index.free) = Q_R(:, pp + i);
        q_R(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q_R(index.node_quasiperiodic1);
        abs_q_R = abs(q_R);
        nodes = nodes_origin + (i - 1) * parameters.lattice.shift; % parameters.lattice.shift;
        patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q_R, ...
            'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
    end

    colorbar
    axis equal
    axis off

else

    
end

end