function fems = applyBoundaryConditions(fems, meshs, parameters)

geometries = fieldnames(meshs);
ngeom = length(geometries);

for gg = 1 : ngeom

    mesh = meshs.(geometries{gg});
    fem = fems.(geometries{gg});
    boundary_conditions = parameters.math.(geometries{gg}).boundary_conditions;

    for i = 1 : length(boundary_conditions)
        switch boundary_conditions{i} 
            case 'Dirichlet'
                fem = applyDirichlet(fem, mesh, parameters);
            case 'Neumann'
                fem = applyNeumann(fem, mesh, parameters);
            case 'quasiperiodic'
                % not applying QBC here.
            case 'Robin'
                fem = applyRobin(fem, mesh, parameters);
            case 'infinite'
                if ismember('Dirichlet', boundary_conditions)
                    fem = applyInfinite(fem, mesh, parameters, 'Dirichlet');
                elseif ismember('Neumann', boundary_conditions)
                    fem = applyInfinite(fem, mesh, parameters, 'Neumann');
                elseif ismember('Robin', boundary_conditions)
                    fem = applyInfinite(fem, mesh, parameters, 'Robin');
                else
                    fem = applyInfinite(fem, mesh, parameters, 'Dirichlet');
                end
            case 'none'
                fem = applyNone(fem, mesh, parameters);
        end
    end

    fem.a_except_quasiperiodic = fem.a;
    fem.b_except_quasiperiodic = fem.b;
    fem.c_except_quasiperiodic = fem.c;
    fem.d_except_quasiperiodic = fem.d;

    fem.a_except_quasiperiodic_0 = fem.a_0;
    fem.b_except_quasiperiodic_0 = fem.b_0;
    fem.c_except_quasiperiodic_0 = fem.c_0;
    fem.d_except_quasiperiodic_0 = fem.d_0;

    fem.a_except_quasiperiodic_end = fem.a_end;
    fem.b_except_quasiperiodic_end = fem.b_end;
    fem.c_except_quasiperiodic_end = fem.c_end;
    fem.d_except_quasiperiodic_end = fem.d_end;

    fems.(geometries{gg}) = fem;

end

end