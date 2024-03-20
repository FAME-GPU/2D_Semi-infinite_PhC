function plotEdgeStates(result, mesh, parameters)

gamma = parameters.gamma;

figure
ev_select = result{parameters.selected_index_wave_vec}.ev_finite(:, parameters.selected_index_eigencurve);
nodes_origin = mesh.geometry1.nodes;
index = mesh.geometry1.index;
elements = mesh.geometry1.elements;
n = mesh.geometry1.n_nodes;
m = mesh.geometry1.n_free_nodes;
if ismember('Dirichlet', parameters.math.geometry1.boundary_conditions)
    ev_select = [zeros(length(index.node_infinite_opposite) - 1, 1); ev_select];
end

q1_all = zeros(n, 1);
q1 = ev_select(1 : m);
q1 = q1 / norm(q1);
q1_all(index.free) = q1;
q1_all(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_all(index.node_quasiperiodic1);
abs_q1_all = abs(q1_all);
nodes = nodes_origin;
patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q1_all, ...
    'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
colorbar
axis equal
axis off

figure
for i = 1 : parameters.num_blocks
    q1_i_all = zeros(n, 1);
    q1_i = ev_select((i - 1) * m + 1 : i * m);
    q1_i_all(index.free) = q1_i;
    q1_i_all(index.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_i_all(index.node_quasiperiodic1);
    abs_q1_i_all = abs(q1_i_all);
    nodes = nodes_origin + (i - 1) * parameters.lattice.shift; % parameters.lattice.shift;
    patch('Faces', elements(1 : 3, :).', 'Vertices', nodes.', 'FaceVertexCData', abs_q1_i_all, ... 
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end
colorbar
axis equal
axis off

end