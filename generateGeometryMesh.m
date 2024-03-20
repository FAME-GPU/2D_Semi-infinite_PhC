function [model, geometry] = generateGeometryMesh(geometry, parameters, flag)
%%

geometries = fieldnames(geometry);
ngeom = length(geometries);

for i = 1 : ngeom
    model.(geometries{i}).model = createpde(1);

    geometryFromEdges(model.(geometries{i}).model, geometry.(geometries{i}).dl);
    geometry.(geometries{i}).model.Geometry = model.(geometries{i}).model.Geometry;
    generateMesh(model.(geometries{i}).model, 'GeometricOrder', 'linear', 'Hmax', parameters.mesh.hmax);  % 'Hvertex', {3:6, 0.1}, 'Hedge', {1:8, 0.1}

    if flag.plot_geometry
        hax = axes(figure);
        pdegplot(geometry.(geometries{i}).dl, 'VertexLabels', 'off', 'EdgeLabels', 'off', 'FaceLabels', 'off', 'FaceAlpha', 1), hold on
        axis equal
        axis off

        switch parameters.lattice.lattice_type_specific
            case 'Square_OME2019'
                arrow([0, 0], parameters.lattice.a1, 'facealpha', 1, 'color', 'r');
                arrow([0, 0], parameters.lattice.a2, 'facealpha', 1, 'color', 'r');
                text(-0.05, -0.05, '0', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                text(0.45 * parameters.lattice.a1(1), parameters.lattice.a1(2) - 0.05, '$\mathbf{a}_1$', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                text(-0.1, parameters.lattice.a2(2) / 2, '$\mathbf{a}_2$', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                tmp_a = (1 - parameters.a_s) / 2;
                plot(hax, [tmp_a, tmp_a + parameters.a_s], [tmp_a, tmp_a], 'Color', 'm', 'LineWidth', 1.5), hold on
                text(0.45 * parameters.lattice.a1(1), parameters.lattice.a1(2) + 0.1, '$a_s$', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'm');
                text(0.45 * parameters.lattice.a1(1), 0.45 * parameters.lattice.a2(2), '$\varepsilon_r$', 'interpret', 'latex', ...
                    'FontSize', 24, 'Color', 'k');
            case 'Hexagonal_PRB2018'
                arrow([0, 0], parameters.lattice.a1, 'facealpha', 1, 'color', 'r');
                arrow([0, 0], parameters.lattice.a2, 'facealpha', 1, 'color', 'r');
                text(-0.05, -0.05, '0', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                text(0.45 * parameters.lattice.a1(1), parameters.lattice.a1(2) - 0.05, '$\mathbf{a}_1$', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                text(0.26 * parameters.lattice.a2(1), parameters.lattice.a2(2) / 2, '$\mathbf{a}_2$', 'interpret', 'latex', ...
                    'FontSize', 18, 'Color', 'r');
                O = (parameters.lattice.a1 + parameters.lattice.a2) / 2;
                plot(hax, [O(1), O(1) + parameters.r1], [O(2), O(2)], 'Color', 'm', 'LineWidth', 1.5), hold on
                text(O(1) + parameters.r1 + 0.05, O(2), '$r_1$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'm');
                plot(hax, [O(1), O(1) + parameters.r2 / 2], [O(2), O(2) - parameters.r2 * sqrt(3) / 2], 'Color', 'k', 'LineWidth', 1.5)
                text(O(1) + parameters.r2 / 2, O(2) - parameters.r2 * sqrt(3) / 2 - 0.05, '$r_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k');
            case 'Hexagonal_PRA2021_2'
                plot(hax, [parameters.lattice.a2(1), parameters.lattice.a2(1) + parameters.lattice.a1(1)], ...
                    [parameters.lattice.a2(2), parameters.lattice.a2(2) + parameters.lattice.a1(2)], ...
                    'LineWidth', 0.5, 'Color', [0, 0.45, 0.74]), hold on
                plot(hax, [2 * parameters.lattice.a2(1), 2 * parameters.lattice.a2(1) + parameters.lattice.a1(1)], ...
                    [2 * parameters.lattice.a2(2), 2 * parameters.lattice.a2(2) + parameters.lattice.a1(2)], ...
                    'LineWidth', 0.5, 'Color', [0, 0.45, 0.74]), hold on
                plot(hax, [3.5 * parameters.lattice.a2(1), 3.5 * parameters.lattice.a2(1) + parameters.lattice.a1(1)], ...
                    [3.5 * parameters.lattice.a2(2), 3.5 * parameters.lattice.a2(2) + parameters.lattice.a1(2)], ...
                    'LineWidth', 0.5, 'Color', [0, 0.45, 0.74]), hold on
                plot(hax, [4.5 * parameters.lattice.a2(1), 4.5 * parameters.lattice.a2(1) + parameters.lattice.a1(1)], ...
                    [4.5 * parameters.lattice.a2(2), 4.5 * parameters.lattice.a2(2) + parameters.lattice.a1(2)], ...
                    'LineWidth', 0.5, 'Color', [0, 0.45, 0.74]), hold on
                plot(hax, [2 * parameters.lattice.a2(1), 2 * parameters.lattice.a2(1) + 0.38 * parameters.lattice.a1(1)], ...
                    [2 * parameters.lattice.a2(2), 2 * parameters.lattice.a2(2)], ...
                    'LineWidth', 1.5, 'Color', 'm'), hold on
                text(2 * parameters.lattice.a2(1) + 0 * parameters.lattice.a1(1), 2 * parameters.lattice.a2(2) + 0.1, '$R$', 'interpret', 'latex', ...
                    'FontSize', 28, 'Color', 'm');
                arrow([0, 0], parameters.lattice.a1, 'facealpha', 1, 'color', 'r');
                arrow([0, 0], parameters.lattice.a2, 'facealpha', 1, 'color', 'r');
                text(-0.1, -0.1, '0', 'interpret', 'latex', ...
                    'FontSize', 28, 'Color', 'r');
                text(0.4 * parameters.lattice.a1(1), parameters.lattice.a1(2) - 0.1, '$\mathbf{a}_1$', 'interpret', 'latex', ...
                    'FontSize', 28, 'Color', 'r');
                text(-0.1 * parameters.lattice.a2(1), parameters.lattice.a2(2) / 2, '$\mathbf{a}_2$', 'interpret', 'latex', ...
                    'FontSize', 28, 'Color', 'r');
                text(0.6, 0.55, '$r_1$', 'interpret', 'latex', 'FontSize', 28, 'Color', 'r');
                text(0.45 + 3 * parameters.lattice.a2(1), 0.0 + 3 * parameters.lattice.a2(2), '$r_2$', 'interpret', 'latex', 'FontSize', 28, 'Color', 'r');
                text(0.3, -0.2, '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
                text(0.3 - 0.1 * parameters.lattice.a2(1), -0.2 - 0.1 * parameters.lattice.a2(2), '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
                text(0.3 - 0.2 * parameters.lattice.a2(1), -0.2 - 0.2 * parameters.lattice.a2(2), '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
                text(0.5 + 5.5 * parameters.lattice.a2(1), 0.2 + 5.5 * parameters.lattice.a2(2), '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
                text(0.5 + 5.6 * parameters.lattice.a2(1), 0.2 + 5.6 * parameters.lattice.a2(2), '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
                text(0.5 + 5.7 * parameters.lattice.a2(1), 0.2 + 5.7 * parameters.lattice.a2(2), '.', 'interpret', 'latex', 'FontSize', 60, 'Color', 'b');
        end

    end

    if flag.plot_mesh
        figure
        pdeplot(model.(geometries{i}).model, 'NodeLabels', 'on', 'ElementLabels', 'off'), hold on;
        axis equal
        axis off

        if strcmp(parameters.lattice.lattice_type_specific, 'sample_semiinfinite')
            a1 = parameters.lattice.a1;
            a2 = parameters.lattice.a2;
            plot([0, a1(1)], [0, a1(2)], 'k-', 'LineWidth', 2);
            plot([a1(1), a1(1) + a2(1)], [a1(2), a1(2) + a2(2)], 'k-', 'LineWidth', 2);
            plot([a1(1) + a2(1), a2(1)], [a1(2) + a2(2), a2(2)], 'k-', 'LineWidth', 2);
            plot([a2(1), 0], [a2(2), 0], 'k-', 'LineWidth', 2);

            a1x = linspace(0, a1(1), 11);
            a1y = linspace(0, a1(2), 11);
            plot(a1x, a1y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
            a2x = linspace(0, a2(1), 8);
            a2y = linspace(0, a2(2), 8);
            plot(a2x, a2y, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
            a2x = 1 + a2x;
            plot(a2x, a2y, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
            a1y = a1y + a2(1);
            a1y = a1y - a2(1) + a2(2);
            a1x = a1x + a2(1);
            plot(a1x, a1y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', '#EDB120', 'MarkerFaceColor', '#EDB120');
        end
        if strcmp(parameters.lattice.lattice_type_specific, 'sample_biinfinite')
            a1 = parameters.lattice.a1;
            a2 = parameters.lattice.a2;
            plot([0, 2 * a1(1)], [0, 2 * a1(2)], 'k-', 'LineWidth', 2);
            plot([2 * a1(1), 2 * a1(1) + a2(1)], [2 * a1(2), 2 * a1(2) + a2(2)], 'k-', 'LineWidth', 2);
            plot([2 * a1(1) + a2(1), a2(1)], [2 * a1(2) + a2(2), a2(2)], 'k-', 'LineWidth', 2);
            plot([a2(1), 0], [a2(2), 0], 'k-', 'LineWidth', 2);
            plot([a1(1), a1(1) + a2(1)], [a1(2), a1(2) + a2(2)], 'k-', 'LineWidth', 2);

            a1x = linspace(0, 2 * a1(1), 21);
            a1y = linspace(0, 2 * a1(2), 21);
            plot(a1x, a1y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
            a2x = linspace(0, a2(1), 11);
            a2y = linspace(0, a2(2), 11);
            plot(a2x, a2y, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
            a2x = a1(1) + a2x;
            plot(a2x, a2y, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
            a2x = a1(1) + a2x;
            plot(a2x, a2y, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
            a1y = a1y + a2(1);
            a1y = a1y - a2(1) + a2(2);
            a1x = a1x + a2(1);
            plot(a1x, a1y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', '#EDB120', 'MarkerFaceColor', '#EDB120');
        end
    end
    geometry.(geometries{i}).model.Mesh = model.(geometries{i}).model.Mesh;

    

        % nodes = model.(geometries{i}).model.Mesh.Nodes;
        % elements = model.(geometries{i}).model.Mesh.Elements;
        % 
        % if 0
        %     % no need to delete nodes now
        % 
        %     % xv = geometry.(geometries{i}).boundary.polygon.vertex(:, 1);
        %     % yv = geometry.(geometries{i}).boundary.polygon.vertex(:, 2);
        %     % [in, on] = inpolygonTol(nodes(:, 1), nodes(:, 2), xv, yv, 1e-8);
        %     % nodes_out = nodes(~in&~on, :);
        % 
        %     [~, out] = isInPolygon(nodes, geometry.(geometries{i}).boundary.polygon.vertex);
        %     [elements, nodes] = RemovePoint(elements', nodes', out', 1);
        %     if 0
        %         figure
        %         trimesh(elements, nodes(:, 1), nodes(:, 2))
        %     end
        %     elements = elements';
        %     nodes = nodes';
        % end
        % 
        % mesh.(geometries{i}).nodes = nodes;
        % mesh.(geometries{i}).elements = elements;
    
    
    
end

end