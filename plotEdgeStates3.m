function plotEdgeStates3(result, mesh, parameters)

gamma = parameters.gamma;

figure
ev_select = result{parameters.selected_index_wave_vec}.ev_finite(:, parameters.selected_index_eigencurve);
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

q1_all = zeros(n2, 1);
q1 = ev_select(parameters.num_blocks * m1 + 1 : parameters.num_blocks * m1 + m2);
q1 = q1 / norm(q1);
q1_all(index2.free) = q1;
q1_all(index2.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_all(index2.node_quasiperiodic1);
abs_q1_all = abs(q1_all);
nodes2 = nodes_origin2;
patch('Faces', elements2(1 : 3, :).', 'Vertices', nodes2.', 'FaceVertexCData', abs_q1_all, ...
    'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
colorbar
axis equal
axis off


figure
% if ismember('infinite', parameters.math.geometry1.boundary_conditions)
%     ev_select = [zeros(length(index1.node_purely_infinite) - 1, 1); ev_select];
% end
% if ismember('infinite', parameters.math.geometry3.boundary_conditions)
%     ev_select = [ev_select; zeros(length(index3.node_purely_infinite) - 1, 1)];
% end
nodes_origin1 = nodes_origin1 - parameters.num_blocks * parameters.lattice.shift1;
ev_select1 = ev_select(1 : parameters.num_blocks * m1);
for i = 1 : parameters.num_blocks
    q1_i_all = zeros(n1, 1);
    q1_i = ev_select1((i - 1) * m1 + 1 : i * m1);
    q1_i_all(index1.free) = q1_i;
    q1_i_all(index1.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q1_i_all(index1.node_quasiperiodic1);
    abs_q1_i_all = abs(q1_i_all);
    nodes1 = nodes_origin1 + (i - 1) * parameters.lattice.shift1; % parameters.lattice.shift;
    patch('Faces', elements1(1 : 3, :).', 'Vertices', nodes1.', 'FaceVertexCData', abs_q1_i_all, ... 
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end

ev_select2 = ev_select(parameters.num_blocks * m1 + 1 : parameters.num_blocks * m1 + m2);
q2_i_all = zeros(n2, 1);
q2_i = ev_select2;
q2_i_all(index2.free) = q2_i;
q2_i_all(index2.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q2_i_all(index2.node_quasiperiodic1);
abs_q2_i_all = abs(q2_i_all);
nodes2 = nodes_origin2; % parameters.lattice.shift;
patch('Faces', elements2(1 : 3, :).', 'Vertices', nodes2.', 'FaceVertexCData', abs_q2_i_all, ... 
    'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on

nodes_origin3 = nodes_origin3 + 1.5 * parameters.lattice.shift1;
ev_select3 = ev_select(parameters.num_blocks * m1 + m2 + 1 : parameters.num_blocks * m1 + m2 + parameters.num_blocks * m3);
for i = 1 : parameters.num_blocks
    q3_i_all = zeros(n3, 1);
    q3_i = ev_select3((i - 1) * m3 + 1 : i * m3);
    q3_i_all(index3.free) = q3_i;
    q3_i_all(index3.node_quasiperiodic2) = gamma(parameters.selected_index_wave_vec) * q3_i_all(index3.node_quasiperiodic1);
    abs_q3_i_all = abs(q3_i_all);
    nodes3 = nodes_origin3 + (i - 1) * parameters.lattice.shift2; % parameters.lattice.shift;
    patch('Faces', elements3(1 : 3, :).', 'Vertices', nodes3.', 'FaceVertexCData', abs_q3_i_all, ... 
        'FaceColor', 'interp', 'EdgeColor', 'interp'), hold on
end
colorbar
axis equal
axis off

end