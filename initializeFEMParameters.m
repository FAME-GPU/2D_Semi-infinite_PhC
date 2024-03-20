function fems = initializeFEMParameters(meshs, parameters)

geometries = fieldnames(meshs);
ngeom = length(geometries);

for gg = 1 : ngeom

    mesh = meshs.(geometries{gg});

    fem.inv_tensor22 = kron(mesh.elements(4, :) ~= 1, parameters.math.inv_tensor22 - eye(2)) + kron(ones(1, mesh.n_elements), eye(2));
    fem.inv_tensor22 = reshape(fem.inv_tensor22, 2, 2, mesh.n_elements);

    fem.tensor11 = kron(mesh.elements(4, :) ~= 1, parameters.math.tensor11 - 1) + kron(ones(1, mesh.n_elements), 1);
    % fem.tensor11 = reshape(fem.tensor11, 1, 1, mesh.n_elements);

    fem.row_indices = reshape(repmat(mesh.elements_simple, 3, 1), [], 1);
    fem.col_indices = reshape(repmat(mesh.elements_simple(:), 1, 3).', 1, []).';

    [fem.rspts, fem.qwgts] = GaussPoints(parameters.fem.Gauss);
    xx = zeros(3, mesh.n_elements);
    yy = zeros(3, mesh.n_elements);
    for m = 1 : 3
        xx(m, :) = mesh.nodes(1, mesh.elements(m, :));
        yy(m, :) = mesh.nodes(2, mesh.elements(m, :));
    end
    [fem.area, fem.a, fem.b, fem.c] = hatGradients(xx, yy);
    fem.d = ones(3, mesh.n_elements);

    fem.n_Gauss_points = size(fem.rspts, 1);
    r = fem.rspts(:, 1);
    s = fem.rspts(:, 2);
    x = mesh.nodes(1, :);
    x = x(mesh.elements(1 : 3, :));
    y = mesh.nodes(2, :);
    y = y(mesh.elements(1 : 3, :));
    [fem.S, fem.detJ] = isoparametricMapVectorize(x, y, r, s, @P1);
    fem.wxarea = (fem.qwgts / 2) .* fem.detJ;
    % fem.SS = zeros(9, length(r), mesh.n_elements);
    % for i = 1 : length(r)
    %     fem.SS(:, :, i) = reshape(fem.S(:, i) * fem.S(:, i)', 9, 1);
    % end
    fem.prod_tensor11_wxarea = reshape(fem.tensor11 .* fem.wxarea, length(r), 1, []);

    fem.a_0 = fem.a;
    fem.b_0 = fem.b;
    fem.c_0 = fem.c;
    fem.d_0 = fem.d;

    fem.a_end = fem.a;
    fem.b_end = fem.b;
    fem.c_end = fem.c;
    fem.d_end = fem.d;

    fems.(geometries{gg}) = fem;

end

end