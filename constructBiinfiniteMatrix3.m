function mtx = constructBiinfiniteMatrix3(fems, meshs, parameters)

geometries = fieldnames(meshs);
ngeom = length(geometries);

if ngeom ~= 3
    error('Number of geometries must be 3.');
end

mesh_L = meshs.(geometries{1});
fem_L = fems.(geometries{1});
mesh_C = meshs.(geometries{2});
fem_C = fems.(geometries{2});
mesh_R = meshs.(geometries{3});
fem_R = fems.(geometries{3});

%% construct left-hand-side matrix
K_C = fem_C.stiffness;
K_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic1) = ...
    K_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic1) + ...
    K_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic2);
K_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_near_quasiperiodic2) = ...
    K_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_near_quasiperiodic2) + ...
    K_C(mesh_C.index.node_quasiperiodic2, mesh_C.index.node_near_quasiperiodic2);
K_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_quasiperiodic1) = ...
    K_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_quasiperiodic1) + ...
    K_C(mesh_C.index.node_quasiperiodic2, mesh_C.index.node_quasiperiodic2);

K_L = fem_L.stiffness;
K_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) = ...
    K_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) + ...
    K_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic2);
K_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) = ...
    K_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) + ...
    K_L(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_near_quasiperiodic2);
K_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) = ...
    K_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) + ...
    K_L(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_quasiperiodic2);

K_LE = fem_L.stiffness_end;
K_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) = ...
    K_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) + ...
    K_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic2);
K_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) = ...
    K_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) + ...
    K_LE(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_near_quasiperiodic2);
K_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) = ...
    K_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) + ...
    K_LE(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_quasiperiodic2);

K_R = fem_R.stiffness;
K_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) = ...
    K_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) + ...
    K_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic2);
K_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) = ...
    K_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) + ...
    K_R(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_near_quasiperiodic2);
K_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) = ...
    K_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) + ...
    K_R(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_quasiperiodic2);

K_RE = fem_R.stiffness_end;
K_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) = ...
    K_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) + ...
    K_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic2);
K_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) = ...
    K_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) + ...
    K_RE(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_near_quasiperiodic2);
K_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) = ...
    K_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) + ...
    K_RE(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_quasiperiodic2);

% K_C(mesh_C.index.node_none{1}, mesh_C.index.node_none{1}) = ...
%     K_C(mesh_C.index.node_none{1}, mesh_C.index.node_none{1}) + ...
%     K_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite);
K_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) = ...
    K_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) + ...
    K_L(mesh_L.index.node_infinite, mesh_L.index.node_infinite);  % maybe add by K_C
K_LE(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) = ...
    K_LE(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) + ...
    K_L(mesh_L.index.node_infinite, mesh_L.index.node_infinite);
K_R(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) = ...
    K_R(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) + ...
    K_R(mesh_R.index.node_infinite, mesh_R.index.node_infinite);
K_RE(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) = ...
    K_RE(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) + ...
    K_R(mesh_R.index.node_infinite, mesh_R.index.node_infinite);

mtx.K_C = K_C(mesh_C.index.free, mesh_C.index.free);
mtx.K_R = K_R(mesh_R.index.free, mesh_R.index.free);
mtx.K_RE = K_RE(mesh_R.index.free_end, mesh_R.index.free_end);
mtx.K_L = K_L(mesh_L.index.free, mesh_L.index.free);
mtx.K_LE = K_LE(mesh_L.index.free_end, mesh_L.index.free_end);  % the left end part has the same DOF as others

E_R = K_R(mesh_R.index.node_purely_near_infinite, mesh_R.index.node_purely_infinite);
% assume E_RB = E_R
mtx.E_RB = sparse(mesh_R.n_free_nodes, mesh_C.n_free_nodes);
mtx.E_R = sparse(mesh_R.n_free_nodes, mesh_R.n_free_nodes);
[n_row_E_R, n_col_E_R] = size(E_R);
mtx.E_RB(1 : n_row_E_R, end - n_col_E_R + 1 : end) = E_R;
mtx.E_R(1 : n_row_E_R, end - n_col_E_R + 1 : end) = E_R;
mtx.E_RE = [mtx.E_R; sparse(length(mesh_R.index.free_end) - length(mesh_R.index.free), length(mesh_R.index.free))];

E_L = K_L(mesh_L.index.node_purely_infinite, mesh_L.index.node_purely_near_infinite);
% assume E_LB = E_L
mtx.E_LB = sparse(mesh_C.n_free_nodes, mesh_L.n_free_nodes);
mtx.E_L = sparse(mesh_L.n_free_nodes, mesh_L.n_free_nodes);
[n_row_E_L, n_col_E_L] = size(E_L);
mtx.E_LB(1 : n_row_E_L, end - n_col_E_L + 1 : end) = E_L;
mtx.E_L(1 : n_row_E_L, end - n_col_E_L + 1 : end) = E_L;
mtx.E_LE = [sparse(mesh_L.n_free_nodes, length(mesh_L.index.free_end) - mesh_L.n_free_nodes), mtx.E_L];

%% construct right-hand-side matrix
M_C = fem_C.mass;
M_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic1) = ...
    M_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic1) + ...
    M_C(mesh_C.index.node_near_quasiperiodic2, mesh_C.index.node_quasiperiodic2);
M_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_near_quasiperiodic2) = ...
    M_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_near_quasiperiodic2) + ...
    M_C(mesh_C.index.node_quasiperiodic2, mesh_C.index.node_near_quasiperiodic2);
M_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_quasiperiodic1) = ...
    M_C(mesh_C.index.node_quasiperiodic1, mesh_C.index.node_quasiperiodic1) + ...
    M_C(mesh_C.index.node_quasiperiodic2, mesh_C.index.node_quasiperiodic2);

M_L = fem_L.mass;
M_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) = ...
    M_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) + ...
    M_L(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic2);
M_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) = ...
    M_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) + ...
    M_L(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_near_quasiperiodic2);
M_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) = ...
    M_L(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) + ...
    M_L(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_quasiperiodic2);

M_LE = fem_L.mass_end;
M_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) = ...
    M_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic1) + ...
    M_LE(mesh_L.index.node_near_quasiperiodic2, mesh_L.index.node_quasiperiodic2);
M_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) = ...
    M_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_near_quasiperiodic2) + ...
    M_LE(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_near_quasiperiodic2);
M_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) = ...
    M_LE(mesh_L.index.node_quasiperiodic1, mesh_L.index.node_quasiperiodic1) + ...
    M_LE(mesh_L.index.node_quasiperiodic2, mesh_L.index.node_quasiperiodic2);

M_R = fem_R.mass;
M_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) = ...
    M_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) + ...
    M_R(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic2);
M_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) = ...
    M_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) + ...
    M_R(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_near_quasiperiodic2);
M_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) = ...
    M_R(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) + ...
    M_R(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_quasiperiodic2);

M_RE = fem_R.mass_end;
M_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) = ...
    M_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic1) + ...
    M_RE(mesh_R.index.node_near_quasiperiodic2, mesh_R.index.node_quasiperiodic2);
M_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) = ...
    M_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_near_quasiperiodic2) + ...
    M_RE(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_near_quasiperiodic2);
M_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) = ...
    M_RE(mesh_R.index.node_quasiperiodic1, mesh_R.index.node_quasiperiodic1) + ...
    M_RE(mesh_R.index.node_quasiperiodic2, mesh_R.index.node_quasiperiodic2);

% M_C(mesh_C.index.node_none{1}, mesh_C.index.node_none{1}) = ...
%     M_C(mesh_C.index.node_none{1}, mesh_C.index.node_none{1}) + ...
%     M_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite);
M_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) = ...
    M_L(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) + ...
    M_L(mesh_L.index.node_infinite, mesh_L.index.node_infinite);  % maybe add by M_C
M_LE(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) = ...
    M_LE(mesh_L.index.node_infinite_opposite, mesh_L.index.node_infinite_opposite) + ...
    M_L(mesh_L.index.node_infinite, mesh_L.index.node_infinite); 
M_R(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) = ...
    M_R(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) + ...
    M_R(mesh_R.index.node_infinite, mesh_R.index.node_infinite);
M_RE(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) = ...
    M_RE(mesh_R.index.node_infinite_opposite, mesh_R.index.node_infinite_opposite) + ...
    M_R(mesh_R.index.node_infinite, mesh_R.index.node_infinite);

mtx.M_C = M_C(mesh_C.index.free, mesh_C.index.free);
mtx.M_R = M_R(mesh_R.index.free, mesh_R.index.free);
mtx.M_RE = M_RE(mesh_R.index.free_end, mesh_R.index.free_end);
mtx.M_L = M_L(mesh_L.index.free, mesh_L.index.free);
mtx.M_LE = M_LE(mesh_L.index.free_end, mesh_L.index.free_end);  % the left end part has the same DOF as others

F_R = M_R(mesh_R.index.node_purely_near_infinite, mesh_R.index.node_purely_infinite);
% assume F_RB = F_R
mtx.F_RB = sparse(mesh_R.n_free_nodes, mesh_C.n_free_nodes);
mtx.F_R = sparse(mesh_R.n_free_nodes, mesh_R.n_free_nodes);
[n_row_F_R, n_col_F_R] = size(F_R);
mtx.F_RB(1 : n_row_F_R, end - n_col_F_R + 1 : end) = F_R;
mtx.F_R(1 : n_row_F_R, end - n_col_F_R + 1 : end) = F_R;
mtx.F_RE = [mtx.F_R; sparse(length(mesh_R.index.free_end) - length(mesh_R.index.free), length(mesh_R.index.free))];

F_L = M_L(mesh_L.index.node_purely_infinite, mesh_L.index.node_purely_near_infinite);
% assume F_LB = F_L
mtx.F_LB = sparse(mesh_C.n_free_nodes, mesh_L.n_free_nodes);
mtx.F_L = sparse(mesh_L.n_free_nodes, mesh_L.n_free_nodes);
[n_row_F_L, n_col_F_L] = size(F_L);
mtx.F_LB(1 : n_row_F_L, end - n_col_F_L + 1 : end) = F_L;
mtx.F_L(1 : n_row_F_L, end - n_col_F_L + 1 : end) = F_L;
mtx.F_LE = [sparse(mesh_L.n_free_nodes, length(mesh_L.index.free_end) - mesh_L.n_free_nodes), mtx.F_L];


%% construct finite matrices
parameters.num_blocks_L = parameters.num_blocks;
parameters.num_blocks_R = parameters.num_blocks;
num_blocks_L = parameters.num_blocks_L;
num_blocks_R = parameters.num_blocks_R;  % 1 for center

mtx.KK_L = kron(speye(num_blocks_L - 1), mtx.K_L);
mtx.KK_R = kron(speye(num_blocks_R - 1), mtx.K_R);
mtx.KK = blkdiag(mtx.K_LE, mtx.KK_L, mtx.K_C, mtx.KK_R, mtx.K_RE);
row_index_E_L = 2 : num_blocks_L - 2;
col_index_E_L = 1 : num_blocks_L - 3;
mtx.EE_L = sparse(row_index_E_L, col_index_E_L, 1, num_blocks_L - 1, num_blocks_L - 1, num_blocks_L - 2);
mtx.EE_L = kron(mtx.EE_L, mtx.E_L);
mtx.EE_L = blkdiag(sparse(size(mtx.K_LE, 1), size(mtx.K_LE, 2)), mtx.EE_L);
mtx.EE_L(size(mtx.K_LE, 1) + 1 : size(mtx.K_LE, 1) + size(mtx.K_L, 1), 1 : size(mtx.K_LE, 2)) = mtx.E_LE;
% mtx.EE_L = blkdiag(mtx.EE_L, sparse(size(mtx.K_C, 1), size(mtx.K_C, 2)));
% mtx.EE_L(end - size(mtx.K_C, 1) + 1 : end, end - size(mtx.K_C, 2) - size(mtx.K_L, 2) + 1 : end - size(mtx.K_C, 2)) = mtx.E_LB;
row_index_E_R = 2 : num_blocks_R - 2;
col_index_E_R = 1 : num_blocks_R - 3;
mtx.EE_R = sparse(row_index_E_R, col_index_E_R, 1, num_blocks_R - 1, num_blocks_R - 1, num_blocks_R - 2);
mtx.EE_R = kron(mtx.EE_R, mtx.E_R);
mtx.EE_R = blkdiag(sparse(size(mtx.K_C, 1), size(mtx.K_C, 2)), mtx.EE_R, sparse(size(mtx.K_RE, 1), size(mtx.K_RE, 2)));
mtx.EE_R(size(mtx.K_C, 1) + 1 : size(mtx.K_C, 1) + size(mtx.K_R, 1), 1 : size(mtx.K_C, 2)) = mtx.E_RB;
mtx.EE_R(end - size(mtx.K_RE, 1) + 1 : end, end - size(mtx.K_RE, 2) - size(mtx.K_R, 2) + 1 : end - size(mtx.K_RE, 2)) = mtx.E_RE;
mtx.EE = blkdiag(mtx.EE_L, mtx.EE_R);
mtx.EE(size(mtx.EE_L, 1) + 1 : size(mtx.EE_L, 1) + size(mtx.K_C, 1), size(mtx.EE_L, 2) - size(mtx.K_L, 2) + 1 : size(mtx.EE_L, 2)) = mtx.E_LB;
mtx.KK = mtx.KK + mtx.EE + mtx.EE';

mtx.MM_L = kron(speye(num_blocks_L - 1), mtx.M_L);
mtx.MM_R = kron(speye(num_blocks_R - 1), mtx.M_R);
mtx.MM = blkdiag(mtx.M_LE, mtx.MM_L, mtx.M_C, mtx.MM_R, mtx.M_RE);
row_index_F_L = 2 : num_blocks_L - 2;
col_index_F_L = 1 : num_blocks_L - 3;
mtx.FF_L = sparse(row_index_F_L, col_index_F_L, 1, num_blocks_L - 1, num_blocks_L - 1, num_blocks_L - 2);
mtx.FF_L = kron(mtx.FF_L, mtx.F_L);
mtx.FF_L = blkdiag(sparse(size(mtx.M_LE, 1), size(mtx.M_LE, 2)), mtx.FF_L);
mtx.FF_L(size(mtx.M_LE, 1) + 1 : size(mtx.M_LE, 1) + size(mtx.M_L, 1), 1 : size(mtx.M_LE, 2)) = mtx.F_LE;
% mtx.FF_L = blkdiag(mtx.FF_L, sparse(size(mtx.M_C, 1), size(mtx.M_C, 2)));
% mtx.FF_L(end - size(mtx.M_C, 1) + 1 : end, end - size(mtx.M_C, 2) - size(mtx.M_L, 2) + 1 : end - size(mtx.M_C, 2)) = mtx.F_LB;
row_index_F_R = 2 : num_blocks_R - 2;
col_index_F_R = 1 : num_blocks_R - 3;
mtx.FF_R = sparse(row_index_F_R, col_index_F_R, 1, num_blocks_R - 1, num_blocks_R - 1, num_blocks_R - 2);
mtx.FF_R = kron(mtx.FF_R, mtx.F_R);
mtx.FF_R = blkdiag(sparse(size(mtx.M_C, 1), size(mtx.M_C, 2)), mtx.FF_R, sparse(size(mtx.M_RE, 1), size(mtx.M_RE, 2)));
mtx.FF_R(size(mtx.M_C, 1) + 1 : size(mtx.M_C, 1) + size(mtx.M_R, 1), 1 : size(mtx.M_C, 2)) = mtx.F_RB;
mtx.FF_R(end - size(mtx.M_RE, 1) + 1 : end, end - size(mtx.M_RE, 2) - size(mtx.M_R, 2) + 1 : end - size(mtx.M_RE, 2)) = mtx.F_RE;
mtx.FF = blkdiag(mtx.FF_L, mtx.FF_R);
mtx.FF(size(mtx.FF_L, 1) + 1 : size(mtx.FF_L, 1) + size(mtx.M_C, 1), size(mtx.FF_L, 2) - size(mtx.M_L, 2) + 1 : size(mtx.FF_L, 2)) = mtx.F_LB;
mtx.MM = mtx.MM + mtx.FF + mtx.FF';

% mtx.KK = blkdiag(mtx.KK(1 : num_blocks_L * size(mtx.K_L, 1), 1 : num_blocks_L * size(mtx.K_L, 2)), ... 
%                  mtx.KK(end - num_blocks_R * size(mtx.K_R, 1) + 1 : end, end - num_blocks_R * size(mtx.K_R, 1) + 1 : end));
% mtx.MM = blkdiag(mtx.MM(1 : num_blocks_L * size(mtx.M_L, 1), 1 : num_blocks_L * size(mtx.M_L, 2)), ... 
%                  mtx.MM(end - num_blocks_R * size(mtx.M_R, 1) + 1 : end, end - num_blocks_R * size(mtx.M_R, 1) + 1 : end));

end