%% main

clear
clear global
clc
close all
warning('off')

nseed = 1;
rng(nseed);

%% Set model parameters
parameters = defaultOptions('Hexagonal_PRA2021', 0.02, 20);
parameters.num_blocks = 3;
switch parameters.lattice.lattice_type_specific
    case 'Square_OME2019'
        parameters.n_eigenvalues = parameters.num_blocks * 5;
        parameters.selected_index_wave_vec = 2;
        parameters.selected_index_eigencurve = 2 * parameters.num_blocks + 1;
    case 'Hexagonal_PRB2018'
        parameters.n_eigenvalues = parameters.num_blocks * 5;
        parameters.selected_index_wave_vec = round(parameters.part_num * 0.5);
        parameters.selected_index_eigencurve = 4 * parameters.num_blocks;
    case 'Hexagonal_PRL2015'
        parameters.n_eigenvalues = parameters.num_blocks * 6;
        parameters.selected_index_wave_vec = 2;
        parameters.selected_index_eigencurve = 2 * parameters.num_blocks + 1;
    case 'Hexagonal_PRA2021'
        parameters.n_eigenvalues = parameters.num_blocks * 10;
        parameters.selected_index_wave_vec = parameters.part_num;
        parameters.selected_index_eigencurve = 2 * parameters.num_blocks + 1;
    case 'Hexagonal_PRA2021_2'
        parameters.n_eigenvalues = parameters.num_blocks * 10;
        parameters.selected_index_wave_vec = parameters.part_num;
        parameters.selected_index_eigencurve = 2 * parameters.num_blocks + 1;
end
parameters.contour_integral.num_points = 12; 
parameters.contour_integral.radius = 2e-2;
parameters.SDA.tol = 1e-12;
parameters.num_propagation = 20;
flag.plot_geometry = 1;
flag.plot_mesh = 0;

parameters = resizeParameters(parameters);

%% Construct geometry
geometry = constructGeometry(parameters.geometry, parameters);

[model, geometry] = generateGeometryMesh(geometry, parameters, flag);

[model, geometry, mesh] = classifyVertices(model, geometry, parameters);

%% Assemble FEM matrices 
% Initialize shape functions
fem = initializeFEMParameters(mesh, parameters);
% apply boundary conditions except for the quasi-periodic boundary condition
fem = applyBoundaryConditions(fem, mesh, parameters);

%% Apply boundary conditions and compute band structures for supercell structures
% (mainly quasi-periodic boundary condition varies)
wave_vec_array = parameters.reciprocal_lattice.wave_vec_array; 
n_wave_vector = parameters.reciprocal_lattice.n_wave_vec;
if strcmp(parameters.lattice.lattice_type_specific, 'Hexagonal_PRA2021')
    gamma = exp(1i * (parameters.lattice.a1' * wave_vec_array));  % no quasi-periodic along a1
elseif strcmp(parameters.lattice.lattice_type_specific, 'Hexagonal_PRA2021_2')
    gamma = exp(1i * (parameters.lattice.a1' * wave_vec_array));  % no quasi-periodic along a1
else
    gamma = exp(1i * (parameters.lattice.a2' * wave_vec_array));  % no quasi-periodic along a1
end
parameters.gamma = gamma;

time_main = tic;
result = cell(n_wave_vector, 1);
for i = 1 : n_wave_vector
    if abs(gamma(i) - 1) < 1e-10
        gamma(i) = exp(1i * 1e-4);
    end
    fem = applyQuasiperiodic(fem, mesh, parameters, gamma(i)); 
    mtx = generateFemMatrices(fem, mesh, parameters);
    % Solve eigenvalue problem
    result{i} = computeEigenpairs(mtx, parameters);
    fprintf('Computing band structure: %2d/%2d. Elapsed time: %.1fs.\n', i, n_wave_vector, result{i}.time_eigs);
end
time_main = toc(time_main);

[frequency, parameters] = plotBandStructures(result, parameters);

if strcmp(parameters.lattice.lattice_type_specific, 'Hexagonal_PRA2021')
    % plot edge states
    plotEdgeStates3(result, mesh, parameters);

    % Employ contour integral to compute selected eigenpairs for infinite structures
    [result_CI, result_select] = computeSelectedEigenpairs3(mtx, result, parameters);

    % Plot results
    % [mtx.P_L, mtx.P_R] = plotFieldPropagation3(result_CI, mtx, mesh, parameters);
    [mtx.P_L, mtx.P_R] = plotFieldPropagation3Split(result_CI, mtx, mesh, parameters);

else
    % plot edge states
    plotEdgeStates(result, mesh, parameters);

    % Employ contour integral to compute selected eigenpairs for infinite structures
    [result_CI, result_select] = computeSelectedEigenpairs(mtx, result, parameters);

    % Plot results
    [propagate_matrix, EmF] = plotFieldPropagation(result_CI, mtx, mesh, parameters);
    [mtx.P_R] = plotFieldPropagationSplit(result_CI, mtx, mesh, parameters);
end

