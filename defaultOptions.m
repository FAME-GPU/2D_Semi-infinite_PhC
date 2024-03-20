function parameters = defaultOptions(lattice_type_specific, hmax, part_num)

parameters.math.is_eigenvalue_problem = 1;
parameters.fem.order = 'P1';
parameters.fem.Gauss = 4;
parameters.part_num = part_num;

switch lattice_type_specific
    case 'Square_OME2019'
        parameters.lattice.lattice_type = 'Square';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        theta = 0;

        parameters.math.equation = 'TE';
        varepsilon_inf = 15.68;
        % epsilon_1 = 1 - (omega_p^2 * (omega + 1i * gamma)) / (omega * ((omega + 1i * gamma)^2 - omega_c^2));
        % epsilon_2 = (omega_p^2 * omega_c) / (omega * ((omega + 1i * gamma)^2 - omega_c^2));
        % epsilon_3 = 1 - omega_p^2 / (omega * (omega + 1i * gamma));
        varepsilon_1 = 0.585;
        varepsilon_2 = 0.190;
        varepsilon_3 = 0.670;
        parameters.math.varepsilon = varepsilon_inf * blkdiag([varepsilon_1, 1i * varepsilon_2; -1i * varepsilon_2, varepsilon_1], varepsilon_3);
        parameters.math.mu = eye(3);  
        % parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'Neumann'};
        parameters.math.geometry1.boundary_conditions = {'infinite', 'quasiperiodic', 'Neumann', 'quasiperiodic'};

        
        a = 1;
        a1 = a * [1; 0];
        a2 = a * [0; 1];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift = -a2;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        % vertices = reciprocal_lattice_vector_b * [-0.5,   0, 0.5;
        %                                              0,   0,   0];
        vertices = 2 * pi * [   0,   0,   0;
                                                  -0.5,   0, 0.5];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);
        
        
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                       a1(1), a1(2); ... 
                                                       a1(1) + a2(1), a1(2) + a2(2); ... 
                                                       a2(1), a2(2)];
        % parameters.geometry.boundary.polygon.vertex(flag_QBC : flag_QBC + 1, 1 + mod(flag_QBC, 2)) = ... 
        %     parameters.geometry.boundary.polygon.vertex(flag_QBC : flag_QBC + 1, 1 + mod(flag_QBC, 2)) - parameters.mesh.erase;
        
        a_s = 0.68 / 2;
        parameters.a_s = a_s * 2;
        parameters.geometry.geometry1.rectangle.rectangle1.vertex = [-a_s * len_a1, -a_s * len_a2; ...
                                                           a_s * len_a1, -a_s * len_a2; ...
                                                           a_s * len_a1, a_s * len_a2; ...
                                                           -a_s * len_a1, a_s * len_a2];
        G = [cos(theta), sin(theta); -sin(theta), cos(theta)];
        parameters.geometry.geometry1.rectangle.rectangle1.vertex = parameters.geometry.geometry1.rectangle.rectangle1.vertex * G';
        parameters.geometry.geometry1.rectangle.rectangle1.vertex = 0.5 * a1.' + 0.5 * a2.' + parameters.geometry.geometry1.rectangle.rectangle1.vertex;
        
        % a_ss = 0.3 / 2;
        % theta2 = pi / 6;
        % parameter.geometry.rectangle.rectangle2.vertex = [-a_ss * len_a1, -a_ss * len_a2; ...
        %                                                   a_ss * len_a1, -a_ss * len_a2; ...
        %                                                   a_ss * len_a1,
        %                                                   a_ss * len_a2; ... 
        %                                                   -a_ss * len_a1, a_ss * len_a2];
        % G = [cos(theta2), sin(theta2); -sin(theta2), cos(theta2)];
        % parameter.geometry.rectangle.rectangle2.vertex = parameter.geometry.rectangle.rectangle2.vertex * G';
        % parameter.geometry.rectangle.rectangle2.vertex = [0.4 * len_a1, 0.4 * len_a2] + parameter.geometry.rectangle.rectangle2.vertex;
    case 'Hexagonal_PRB2018'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        theta = 0;

        parameters.math.equation = 'TM';
        varepsilon_r = 13;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        kappa = 0.4;
        parameters.math.mu = [mu_r, 1i * kappa; -1i * kappa, mu_r];
        parameters.math.mu = blkdiag(parameters.math.mu, 1);
        % parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'Dirichlet'};
        parameters.math.geometry1.boundary_conditions = {'infinite', 'quasiperiodic', 'Dirichlet', 'quasiperiodic'};
        
        a = 1;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift = -a2;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        % vertices = reciprocal_lattice_vector_b * [-0.5,   0, 0.5;
        %                                              0,   0,   0];
        vertices = 2 * pi * [   0,   0,   0;
                                                  -0.5,   0, 0.5];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);
        
        
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                       a1(1), a1(2); ... 
                                                       a1(1) + a2(1), a1(2) + a2(2); ... 
                                                       a2(1), a2(2)];
        % parameters.geometry.boundary.polygon.vertex(flag_QBC : flag_QBC + 1, 1 + mod(flag_QBC, 2)) = ... 
        %     parameters.geometry.boundary.polygon.vertex(flag_QBC : flag_QBC + 1, 1 + mod(flag_QBC, 2)) - parameters.mesh.erase;
        
        r1 = 0.3075 * a;
        r2 = 0.0615 * a;
        parameters.r1 = r1;
        parameters.r2 = r2;
        G1 = [cos(theta), sin(theta); -sin(theta), cos(theta)];
        parameters.geometry.geometry1.polygon.triangle1.vertex = [r1, 0; r2 / 2, r2 * sqrt(3) / 2; 0, 0] * G1.';
        parameters.geometry.geometry1.polygon.triangle2.vertex = [r1, 0; r2 / 2, -r2 * sqrt(3) / 2; 0, 0] * G1.';
        G2 = [-1 / 2, sqrt(3) / 2; -sqrt(3) / 2, -1 / 2];
        parameters.geometry.geometry1.polygon.triangle3.vertex = parameters.geometry.geometry1.polygon.triangle1.vertex * G2.';
        parameters.geometry.geometry1.polygon.triangle4.vertex = parameters.geometry.geometry1.polygon.triangle2.vertex * G2.';
        parameters.geometry.geometry1.polygon.triangle5.vertex = parameters.geometry.geometry1.polygon.triangle3.vertex * G2.';
        parameters.geometry.geometry1.polygon.triangle6.vertex = parameters.geometry.geometry1.polygon.triangle4.vertex * G2.';
        parameters.geometry.geometry1.polygon.triangle1.vertex = parameters.geometry.geometry1.polygon.triangle1.vertex + 0.5 * a1.' + 0.5 * a2.';
        parameters.geometry.geometry1.polygon.triangle2.vertex = parameters.geometry.geometry1.polygon.triangle2.vertex + 0.5 * a1.' + 0.5 * a2.';
        parameters.geometry.geometry1.polygon.triangle3.vertex = parameters.geometry.geometry1.polygon.triangle3.vertex + 0.5 * a1.' + 0.5 * a2.';
        parameters.geometry.geometry1.polygon.triangle4.vertex = parameters.geometry.geometry1.polygon.triangle4.vertex + 0.5 * a1.' + 0.5 * a2.';
        parameters.geometry.geometry1.polygon.triangle5.vertex = parameters.geometry.geometry1.polygon.triangle5.vertex + 0.5 * a1.' + 0.5 * a2.';
        parameters.geometry.geometry1.polygon.triangle6.vertex = parameters.geometry.geometry1.polygon.triangle6.vertex + 0.5 * a1.' + 0.5 * a2.';

    case 'Hexagonal_PRL2015'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 11.7;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        parameters.math.mu = mu_r * eye(3);
        % determine the boundary conditions for the finite system,
        % 'Dirichlet' will be replaced by infinite
        parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'infinite'}; 
        parameters.math.geometry2.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'none'};

        % parameters.math.geometry1.boundary_condition.Dirichlet = 0;
        
        a = 1;
        % a = a * sqrt(3) / 3;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift1 = a1;
        parameters.lattice.shift2 = -a1;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        vertices = reciprocal_lattice_vector_b * [   0,   0,   0;
                                                  -0.5,   0, 0.5];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);

        % R1 = 1 / 2.9;
        % R2 = 1 / 3.125;
        % r = 0.24 / 2;

        R1 = 1 / 3;
        R2 = 1 / 3;
        r = R1 / 3;
                
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry1.circle.circle1.center = (a1 + a2) * R1;
        parameters.geometry.geometry1.circle.circle1.radius = r;
        parameters.geometry.geometry1.circle.circle2.center = (a1 + a2) * (1 - R1);
        parameters.geometry.geometry1.circle.circle2.radius = r;
        parameters.geometry.geometry2.boundary.polygon.n = 4;
        parameters.geometry.geometry2.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry2.circle.circle1.center = (a1 + a2) * R2;
        parameters.geometry.geometry2.circle.circle1.radius = r;
        parameters.geometry.geometry2.circle.circle2.center = (a1 + a2) * (1 - R2);
        parameters.geometry.geometry2.circle.circle2.radius = r;

    case 'Hexagonal_PRA2021'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 9.6;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        parameters.math.mu = mu_r * eye(3);
        % determine the boundary conditions for the finite system,
        % 'Dirichlet' will be replaced by infinite
        parameters.math.geometry1.boundary_conditions = {'infinite', 'quasiperiodic', 'none', 'quasiperiodic'};
        parameters.math.geometry2.boundary_conditions = {'none', 'quasiperiodic', 'none', 'quasiperiodic'};
        parameters.math.geometry3.boundary_conditions = {'none', 'quasiperiodic', 'infinite', 'quasiperiodic'};
        % parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'none'};
        % parameters.math.geometry2.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'none'};
        % parameters.math.geometry3.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'infinite'};

        % parameters.math.geometry1.boundary_condition.Dirichlet = 0;
        
        a = 1;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift1 = a2;
        parameters.lattice.shift2 = a2;
        parameters.mesh.hmax = hmax * len_a1; 
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        vertices = reciprocal_lattice_vector_b * [-0.2, 0, 0.2;
                                                     0, 0,   0];
        % vertices = reciprocal_lattice_vector_b * [   0,   0,   0;
        %                                           -0.4,   0, 0.4];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);

        R = 0.38;
        r = 0.1;
        rs = 0.28;
        
        
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry1.circle.circle1.center = R * a1;
        parameters.geometry.geometry1.circle.circle1.radius = r;
        parameters.geometry.geometry1.circle.circle2.center = (1 - R) * a1;
        parameters.geometry.geometry1.circle.circle2.radius = r;
        parameters.geometry.geometry1.circle.circle3.center = R * a1 + a2;
        parameters.geometry.geometry1.circle.circle3.radius = r;
        parameters.geometry.geometry1.circle.circle4.center = (1 - R) * a1 + a2;
        parameters.geometry.geometry1.circle.circle4.radius = r;
        parameters.geometry.geometry1.circle.circle5.center = R * a2;
        parameters.geometry.geometry1.circle.circle5.radius = r;
        parameters.geometry.geometry1.circle.circle6.center = (1 - R) * a2;
        parameters.geometry.geometry1.circle.circle6.radius = r;
        parameters.geometry.geometry1.circle.circle7.center = R * a2 + a1;
        parameters.geometry.geometry1.circle.circle7.radius = r;
        parameters.geometry.geometry1.circle.circle8.center = (1 - R) * a2 + a1;
        parameters.geometry.geometry1.circle.circle8.radius = r;
        parameters.geometry.geometry1.circle.circle9.center = R * (a2 - a1) + a1;
        parameters.geometry.geometry1.circle.circle9.radius = r;
        parameters.geometry.geometry1.circle.circle10.center = (1 - R) * (a2 - a1) + a1;
        parameters.geometry.geometry1.circle.circle10.radius = r;
        parameters.geometry.geometry1.polygon.polygon1.n = 4;
        parameters.geometry.geometry1.polygon.polygon1.vertex = repmat(-a2', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon1.deleted_by = {'circle1 ', 'circle2 '}; % preserve blank to distinguish from circle10 ...
        parameters.geometry.geometry1.polygon.polygon2.n = 4;
        parameters.geometry.geometry1.polygon.polygon2.vertex = repmat(a2', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon2.deleted_by = {'circle3 ', 'circle4 '};
        parameters.geometry.geometry1.polygon.polygon3.n = 4;
        parameters.geometry.geometry1.polygon.polygon3.vertex = repmat(-a1', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon3.deleted_by = {'circle5 ', 'circle6 '};
        parameters.geometry.geometry1.polygon.polygon4.n = 4;
        parameters.geometry.geometry1.polygon.polygon4.vertex = repmat(a1', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon4.deleted_by = {'circle7 ', 'circle8 '};
        % parameters.geometry.geometry1.formula = ['(circle1 - polygon1 + circle2 - polygon1 + circle3 - polygon2 + circle3 - polygon2', ...
        %                                         '+ circle3 - polygon2 + circle3 - polygon2 + circle3 - polygon2 + circle3 - polygon2', ...
        %                                         '+ circle3 - polygon2 + circle3 - polygon2 + circle3 - polygon2 + circle3 - polygon2'];
        % description = zeros(max_nrow, 18);
        % description(:, 1) = [1; circle1.center; cicle1.radius]
        % 
        % [1; geom.center(:); geom.radius; zeros(max_nrow - 4, 1)];

        parameters.geometry.geometry2.boundary.polygon.n = 4;
        parameters.geometry.geometry2.boundary.polygon.vertex = [0, 0; ... 
                                                       a1(1), a1(2); ... 
                                                       a1(1) + 1.5 * a2(1), a1(2) + 1.5 * a2(2); ... 
                                                       1.5 * a2(1), 1.5 * a2(2)];
        parameters.geometry.geometry2.circle.circle1.center = R * a1;
        parameters.geometry.geometry2.circle.circle1.radius = r;
        parameters.geometry.geometry2.circle.circle2.center = (1 - R) * a1;
        parameters.geometry.geometry2.circle.circle2.radius = r;
        parameters.geometry.geometry2.circle.circle3.center = R * a2;
        parameters.geometry.geometry2.circle.circle3.radius = r;
        parameters.geometry.geometry2.circle.circle4.center = (1 - R) * a2;
        parameters.geometry.geometry2.circle.circle4.radius = r;
        parameters.geometry.geometry2.circle.circle5.center = R * a2 + a1;
        parameters.geometry.geometry2.circle.circle5.radius = r;
        parameters.geometry.geometry2.circle.circle6.center = (1 - R) * a2 + a1;
        parameters.geometry.geometry2.circle.circle6.radius = r;
        parameters.geometry.geometry2.circle.circle7.center = R * (a2 - a1) + a1;
        parameters.geometry.geometry2.circle.circle7.radius = r;
        parameters.geometry.geometry2.circle.circle8.center = (1 - R) * (a2 - a1) + a1;
        parameters.geometry.geometry2.circle.circle8.radius = r;
        parameters.geometry.geometry2.circle.circle9.center = 0.5 * a1 + a2;
        parameters.geometry.geometry2.circle.circle9.radius = rs;
        parameters.geometry.geometry2.polygon.polygon1.n = 4;
        parameters.geometry.geometry2.polygon.polygon1.vertex = repmat(-1.5 * a2', 4, 1) + ... 
            parameters.geometry.geometry2.boundary.polygon.vertex;
        parameters.geometry.geometry2.polygon.polygon1.deleted_by = {'circle1 ', 'circle2 '};
        parameters.geometry.geometry2.polygon.polygon2.n = 4;
        parameters.geometry.geometry2.polygon.polygon2.vertex = repmat(-a1', 4, 1) + ... 
            parameters.geometry.geometry2.boundary.polygon.vertex;
        parameters.geometry.geometry2.polygon.polygon2.deleted_by = {'circle3 ', 'circle4 '};
        parameters.geometry.geometry2.polygon.polygon3.n = 4;
        parameters.geometry.geometry2.polygon.polygon3.vertex = repmat(a1', 4, 1) + ... 
            parameters.geometry.geometry2.boundary.polygon.vertex;
        parameters.geometry.geometry2.polygon.polygon3.deleted_by = {'circle5 ', 'circle6 '};

        parameters.geometry.geometry3.boundary.polygon.n = 4;
        parameters.geometry.geometry3.boundary.polygon.vertex = [0, 0; ... 
                                                       a1(1), a1(2); ... 
                                                       a1(1) + a2(1), a1(2) + a2(2); ... 
                                                       a2(1), a2(2)];
        parameters.geometry.geometry3.circle.circle1.center = (a1 + a2) / 2;
        parameters.geometry.geometry3.circle.circle1.radius = rs;

    case 'Hexagonal_PRA2021_2'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 9.6;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        parameters.math.mu = mu_r * eye(3);
        % determine the boundary conditions for the finite system,
        % 'Dirichlet' will be replaced by infinite
        parameters.math.geometry1.boundary_conditions = {'infinite', 'quasiperiodic', 'none', 'quasiperiodic'};
        parameters.math.geometry2.boundary_conditions = {'none', 'quasiperiodic', 'none', 'quasiperiodic'};
        parameters.math.geometry3.boundary_conditions = {'none', 'quasiperiodic', 'infinite', 'quasiperiodic'};
        % parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'none'};
        % parameters.math.geometry2.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'none'};
        % parameters.math.geometry3.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'infinite'};

        % parameters.math.geometry1.boundary_condition.Dirichlet = 0;
        
        a = 1;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift1 = a2;
        parameters.lattice.shift2 = a2;
        parameters.mesh.hmax = hmax * len_a1; 
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        vertices = reciprocal_lattice_vector_b * [-0.2, 0, 0.2;
                                                     0, 0,   0];
        % vertices = reciprocal_lattice_vector_b * [   0,   0,   0;
        %                                           -0.4,   0, 0.4];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);

        R = 0.38;
        r = 0.1;
        rs = 0.28;
        
        
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + 5.5 * a2(1), a1(2) + 5.5 * a2(2); ... 
                                                                 5.5 * a2(1), 5.5 * a2(2)];

        parameters.geometry.geometry1.circle.circle1.center = R * a1;
        parameters.geometry.geometry1.circle.circle1.radius = r;
        parameters.geometry.geometry1.circle.circle2.center = (1 - R) * a1;
        parameters.geometry.geometry1.circle.circle2.radius = r;
        parameters.geometry.geometry1.circle.circle3.center = R * a1 + a2;
        parameters.geometry.geometry1.circle.circle3.radius = r;
        parameters.geometry.geometry1.circle.circle4.center = (1 - R) * a1 + a2;
        parameters.geometry.geometry1.circle.circle4.radius = r;
        parameters.geometry.geometry1.circle.circle5.center = R * a2;
        parameters.geometry.geometry1.circle.circle5.radius = r;
        parameters.geometry.geometry1.circle.circle6.center = (1 - R) * a2;
        parameters.geometry.geometry1.circle.circle6.radius = r;
        parameters.geometry.geometry1.circle.circle7.center = R * a2 + a1;
        parameters.geometry.geometry1.circle.circle7.radius = r;
        parameters.geometry.geometry1.circle.circle8.center = (1 - R) * a2 + a1;
        parameters.geometry.geometry1.circle.circle8.radius = r;
        parameters.geometry.geometry1.circle.circle9.center = R * (a2 - a1) + a1;
        parameters.geometry.geometry1.circle.circle9.radius = r;
        parameters.geometry.geometry1.circle.circle10.center = (1 - R) * (a2 - a1) + a1;
        parameters.geometry.geometry1.circle.circle10.radius = r;

        parameters.geometry.geometry1.circle.circle11.center = R * a1 + a2 + a2;
        parameters.geometry.geometry1.circle.circle11.radius = r;
        parameters.geometry.geometry1.circle.circle12.center = (1 - R) * a1 + a2 + a2;
        parameters.geometry.geometry1.circle.circle12.radius = r;
        parameters.geometry.geometry1.circle.circle13.center = R * a2 + a2;
        parameters.geometry.geometry1.circle.circle13.radius = r;
        parameters.geometry.geometry1.circle.circle14.center = (1 - R) * a2 + a2;
        parameters.geometry.geometry1.circle.circle14.radius = r;
        parameters.geometry.geometry1.circle.circle15.center = R * a2 + a1 + a2;
        parameters.geometry.geometry1.circle.circle15.radius = r;
        parameters.geometry.geometry1.circle.circle16.center = (1 - R) * a2 + a1 + a2;
        parameters.geometry.geometry1.circle.circle16.radius = r;
        parameters.geometry.geometry1.circle.circle17.center = R * (a2 - a1) + a1 + a2;
        parameters.geometry.geometry1.circle.circle17.radius = r;
        parameters.geometry.geometry1.circle.circle18.center = (1 - R) * (a2 - a1) + a1 + a2;
        parameters.geometry.geometry1.circle.circle18.radius = r;

        parameters.geometry.geometry1.circle.circle19.center = R * a2 + 2 * a2;
        parameters.geometry.geometry1.circle.circle19.radius = r;
        parameters.geometry.geometry1.circle.circle20.center = (1 - R) * a2 + 2 * a2;
        parameters.geometry.geometry1.circle.circle20.radius = r;
        parameters.geometry.geometry1.circle.circle21.center = R * a2 + a1 + 2 * a2;
        parameters.geometry.geometry1.circle.circle21.radius = r;
        parameters.geometry.geometry1.circle.circle22.center = (1 - R) * a2 + a1 + 2 * a2;
        parameters.geometry.geometry1.circle.circle22.radius = r;
        parameters.geometry.geometry1.circle.circle23.center = R * (a2 - a1) + a1 + 2 * a2;
        parameters.geometry.geometry1.circle.circle23.radius = r;
        parameters.geometry.geometry1.circle.circle24.center = (1 - R) * (a2 - a1) + a1 + 2 * a2;
        parameters.geometry.geometry1.circle.circle24.radius = r;
        parameters.geometry.geometry1.circle.circle25.center = 0.5 * a1 + a2 + 2 * a2;
        parameters.geometry.geometry1.circle.circle25.radius = rs;

        parameters.geometry.geometry1.circle.circle26.center = (a1 + a2) / 2 + 3.5 * a2;
        parameters.geometry.geometry1.circle.circle26.radius = rs;

        parameters.geometry.geometry1.circle.circle27.center = (a1 + a2) / 2 + 4.5 * a2;
        parameters.geometry.geometry1.circle.circle27.radius = rs;

        parameters.geometry.geometry1.polygon.polygon1.n = 4;
        parameters.geometry.geometry1.polygon.polygon1.vertex = repmat(-5.5 * a2', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon1.deleted_by = {'circle1 ', 'circle2 '};
        parameters.geometry.geometry1.polygon.polygon2.n = 4;
        parameters.geometry.geometry1.polygon.polygon2.vertex = repmat(-a1', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon2.deleted_by = {'circle5 ', 'circle6 ', 'circle13 ', 'circle14 ', 'circle19 ', 'circle20 '};
        parameters.geometry.geometry1.polygon.polygon3.n = 4;
        parameters.geometry.geometry1.polygon.polygon3.vertex = repmat(a1', 4, 1) + ... 
            parameters.geometry.geometry1.boundary.polygon.vertex;
        parameters.geometry.geometry1.polygon.polygon3.deleted_by = {'circle7 ', 'circle8 ', 'circle15 ', 'circle16 ', 'circle21 ', 'circle22 '};

    case 'Hexagonal_PRL2015_2'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 11.7;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        parameters.math.mu = mu_r * eye(3);
        % determine the boundary conditions for the finite system,
        % 'Dirichlet' will be replaced by infinite
        parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'infinite'}; 
        parameters.math.geometry2.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'none'};

        % parameters.math.geometry1.boundary_condition.Dirichlet = 0;
        
        a = 1;
        % a = a * sqrt(3) / 3;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift1 = a1;
        parameters.lattice.shift2 = -a1;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2])';
        vertices = reciprocal_lattice_vector_b * [   0,   0,   0;
                                                  -0.5,   0, 0.5];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);

        % R1 = 1 / 2.9;
        % R2 = 1 / 3.125;
        % r = 0.24 / 2;

        R1 = 1 / 3;
        R2 = 1 / 3;
        r = R1 / 3;
                
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry1.circle.circle1.center = (a1 + a2) * R1;
        parameters.geometry.geometry1.circle.circle1.radius = r;
        parameters.geometry.geometry1.circle.circle2.center = (a1 + a2) * (1 - R1);
        parameters.geometry.geometry1.circle.circle2.radius = r;
        parameters.geometry.geometry2.boundary.polygon.n = 4;
        parameters.geometry.geometry2.boundary.polygon.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry2.circle.circle1.center = (a1 + a2) * R2;
        parameters.geometry.geometry2.circle.circle1.radius = r;
        parameters.geometry.geometry2.circle.circle2.center = (a1 + a2) * (1 - R2);
        parameters.geometry.geometry2.circle.circle2.radius = r;

    case 'sample_semiinfinite'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 13;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        kappa = 0.4;
        parameters.math.mu = [mu_r, 1i * kappa; -1i * kappa, mu_r];
        parameters.math.mu = blkdiag(parameters.math.mu, 1);
        parameters.math.boundary_condition.self = {'quasiperiodic', 'infinite', 'quasiperiodic', 'Neumann'};
        parameters.math.boundary_condition.Dirichlet = 0;
        
        a = 1;
        a1 = a * [1; 0];
        a2 = a * [1/2; 1/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift = a1;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        % reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2]);
        vertices = 2 * pi * [-0.5, 0, 0.5;
                                0, 0,   0];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);
        
        
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                       a1(1), a1(2); ... 
                                                       a1(1) + a2(1), a1(2) + a2(2); ... 
                                                       a2(1), a2(2)];
        parameters.geometry.geometry1.circle.circle1.center = [(a1(1) + a2(1)) / 2; (a1(2) + a2(2)) / 2];
        parameters.geometry.geometry1.circle.circle1.radius = 0.15 * a;

        % parameters.geometry.geometry1.polygon.polygon1.n = 4;
        % parameters.geometry.geometry1.polygon.polygon1.vertex = [0, 0; ... 
        %                                                a1(1), a1(2); ... 
        %                                                a1(1) + a2(1), a1(2) + a2(2); ... 
        %                                                a2(1), a2(2)] / 2;
        % parameters.geometry.geometry1.polygon.polygon1.delete = 1;

    case 'sample_biinfinite'
        parameters.lattice.lattice_type = 'Hexagonal';
        parameters.lattice.lattice_type_specific = lattice_type_specific;
        parameters.math.equation = 'TM';
        varepsilon_r = 11.7;
        parameters.math.varepsilon = varepsilon_r * eye(3);
        mu_r = 1;
        parameters.math.mu = mu_r * eye(3);
        % determine the boundary conditions for the finite system,
        % 'Dirichlet' will be replaced by infinite
        parameters.math.geometry1.boundary_conditions = {'quasiperiodic', 'none', 'quasiperiodic', 'infinite'}; 
        parameters.math.geometry2.boundary_conditions = {'quasiperiodic', 'infinite', 'quasiperiodic', 'none'};

        % parameters.math.geometry1.boundary_condition.Dirichlet = 0;
        
        a = 1;
        a = a * sqrt(3) / 3;
        a1 = a * [1; 0];
        a2 = a * [1/2; sqrt(3)/2];
        len_a1 = norm(a1);
        len_a2 = norm(a2);
        parameters.lattice.a = a;
        parameters.lattice.a1 = a1;
        parameters.lattice.a2 = a2;
        parameters.lattice.len_a1 = len_a1;
        parameters.lattice.len_a2 = len_a2;
        parameters.lattice.shift1 = a1;
        parameters.lattice.shift2 = -a1;
        parameters.mesh.hmax = hmax * len_a1;
        parameters.mesh.erase = 0.5 * parameters.mesh.hmax;
        parameters.brillouin_zone_point = 'GMXG';
        % reciprocal_lattice_vector_b = 2 * pi * inv([a1, a2]);
        vertices = 2 * pi * [-0.5, 0, 0.5;
                                0, 0,   0];
        parameters.reciprocal_lattice.vertices = vertices;
        wave_vec_array = vertices(:, 1);
        for i = 1 : size(vertices, 2) - 1
            tmp1 = linspace(vertices(1, i), vertices(1, i + 1), part_num);
            tmp2 = linspace(vertices(2, i), vertices(2, i + 1), part_num);
            wave_vec_array = [wave_vec_array, [tmp1(2 : end); tmp2(2 : end)]];
        end
        parameters.reciprocal_lattice.wave_vec_array = wave_vec_array;
        parameters.reciprocal_lattice.n_wave_vec = size(wave_vec_array, 2);

        r = 0.24 / 2;
        R1 = 1 / 2.9;
        R2 = 1 / 3.125;
                
        parameters.geometry.geometry1.boundary.polygon.n = 4;
        parameters.geometry.geometry1.boundary.polygon.vertex = [0, 0; ... 
                                                                 2 * a1(1), 2 * a1(2); ... 
                                                                 2 * a1(1) + a2(1), 2 * a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
        parameters.geometry.geometry1.circle.circle1.center = a1 + (a1 + a2) / 2;
        parameters.geometry.geometry1.circle.circle1.radius = r;
        parameters.geometry.geometry1.polygon.polygon1.n = 4;
        parameters.geometry.geometry1.polygon.polygon1.vertex = [(a1(1) + a2(1)) / 2 - r, (a1(2) + a2(2)) / 2 - r; ... 
                                                               (a1(1) + a2(1)) / 2 + r, (a1(2) + a2(2)) / 2 - r; ... 
                                                               (a1(1) + a2(1)) / 2 + r, (a1(2) + a2(2)) / 2 + r; ... 
                                                               (a1(1) + a2(1)) / 2 - r, (a1(2) + a2(2)) / 2 + r];
        parameters.geometry.geometry1.polygon.polygon2.n = 4;
        parameters.geometry.geometry1.polygon.polygon2.vertex = [0, 0; ... 
                                                                 a1(1), a1(2); ... 
                                                                 a1(1) + a2(1), a1(2) + a2(2); ... 
                                                                 a2(1), a2(2)];
end

end