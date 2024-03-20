function mtx = generateFemMatrices(fems, meshs, parameters)


geometries = fieldnames(meshs);
ngeom = length(geometries);

for gg = 1 : ngeom

    mesh = meshs.(geometries{gg});
    fem = fems.(geometries{gg});
    boundary_conditions = parameters.math.(geometries{gg}).boundary_conditions;

    fem.bc = zeros(3, 2 * mesh.n_elements);
    fem.bc(:, 1 : 2 : 2 * mesh.n_elements - 1) = fem.b;
    fem.bc(:, 2 : 2 : 2 * mesh.n_elements) = fem.c;
    fem.bc = reshape(fem.bc, 3, 2, mesh.n_elements);  % tensor

    fem.bc_0 = zeros(3, 2 * mesh.n_elements);
    fem.bc_0(:, 1 : 2 : 2 * mesh.n_elements - 1) = fem.b_0;
    fem.bc_0(:, 2 : 2 : 2 * mesh.n_elements) = fem.c_0;
    fem.bc_0 = reshape(fem.bc_0, 3, 2, mesh.n_elements);  % tensor

    fem.bc_end = zeros(3, 2 * mesh.n_elements);
    fem.bc_end(:, 1 : 2 : 2 * mesh.n_elements - 1) = fem.b_end;
    fem.bc_end(:, 2 : 2 : 2 * mesh.n_elements) = fem.c_end;
    fem.bc_end = reshape(fem.bc_end, 3, 2, mesh.n_elements);  % tensor

    % if ismember('Dirichlet', boundary_conditions)
    %     fem.bc_infinite_opposite = zeros(3, 2 * mesh.n_elements);
    %     fem.bc_infinite_opposite(:, 1 : 2 : 2 * mesh.n_elements - 1) = fem.b_infinite_opposite;
    %     fem.bc_infinite_opposite(:, 2 : 2 : 2 * mesh.n_elements) = fem.c_infinite_opposite;
    %     fem.bc_infinite_opposite = reshape(fem.bc_infinite_opposite, 3, 2, mesh.n_elements);  % tensor
    % 
    %     fem.bc_infinite = zeros(3, 2 * mesh.n_elements);
    %     fem.bc_infinite(:, 1 : 2 : 2 * mesh.n_elements - 1) = fem.b_infinite;
    %     fem.bc_infinite(:, 2 : 2 : 2 * mesh.n_elements) = fem.c_infinite;
    %     fem.bc_infinite = reshape(fem.bc_infinite, 3, 2, mesh.n_elements);  % tensor
    % end

    switch parameters.fem.order
        case 'P1'
            fem = generateStiffnessMatrix(fem, mesh, boundary_conditions);
            fem = generateMassMatrix(fem, mesh, boundary_conditions);

        case 2
            fem = 2;
    end

    fems.(geometries{gg}) = fem;

end

if ngeom == 1
    mtx = constructSemiinfiniteMatrix(fems, meshs, parameters);
elseif ngeom == 2
    mtx = constructBiinfiniteMatrix2(fems, meshs, parameters);
elseif ngeom == 3
    mtx = constructBiinfiniteMatrix3(fems, meshs, parameters);
end

end