function geometry = constructGeometry(geometry, parameters)

geometries = fieldnames(geometry);
switch parameters.lattice.lattice_type
    case 'Square'
        max_nrow = 10;
    case 'Hexagonal'
        max_nrow = 8;
        switch parameters.lattice.lattice_type_specific
            case 'Hexagonal_PRA2021'
                max_nrow = 4;
        end
end
max_nrow = 10;
ngeom = length(geometries);

for k = 1 : ngeom
    nshapes = 0;
    description = [];
    formula = '';
    shapes = fieldnames(geometry.(geometries{k}));
    for i = 1 : length(shapes)
        shapes_i = fieldnames(geometry.(geometries{k}).(shapes{i})); 
        for j = 1 : length(shapes_i)
            geom = geometry.(geometries{k}).(shapes{i}).(shapes_i{j});
            switch shapes{i}
                case 'circle'
                    eval_string = [shapes_i{j} ' = [1; geom.center(:); geom.radius; zeros(max_nrow - 4, 1)];'];
                case 'polygon'
                    eval_string = [shapes_i{j} ' = [2; size(geom.vertex, 1); geom.vertex(:); zeros(max_nrow - (2 * size(geom.vertex, 1) + 2), 1)];'];
                case 'rectangle'
                    eval_string = [shapes_i{j} ' = [3; 4; geom.vertex(:); zeros(max_nrow - 10, 1)];'];
                case 'ellipse'
                    eval_string = [shapes_i{j} ' = [4; geom.center(:); geom.semiaxes(:); geom.angle; zeros(max_nrow - 6, 1)];'];
                case 'boundary'
                    eval_string = [shapes_i{j} ' = [2; geom.n; geom.vertex(:); zeros(max_nrow - (2 * geom.n + 2), 1)];'];
            end
            eval(eval_string);
            
            eval(['description = [description, ' shapes_i{j} '];']);
            if i == 1 && j == 1
                names = char(shapes_i{j});
            else
                names = char(names, shapes_i{j});
            end
            if isfield(geom, 'deleted_by')
                % formula = append(formula, ' - ');
                if ~isempty(geom.deleted_by)
                    for m = 1 : length(geom.deleted_by)
                        formula = insertAfter(formula, geom.deleted_by{m}, ['- ', shapes_i{j}]);
                    end
                end
            else
                formula = append(formula, ' + ');
                formula = append(formula, shapes_i{j});
            end
            
            nshapes = nshapes + 1;
        end
    end

    if isfield(geometry.(geometries{k}), 'boundary')
        if 0
            geometry.(geometries{k}).n = nshapes - 1;
            geometry.(geometries{k}).description = description(:, 2 : end);
            geometry.(geometries{k}).names = names(2 : end, :)';
            formula = ['(', formula(14 : end)];
            geometry.(geometries{k}).formula = append(formula, ')');
            [geometry.(geometries{k}).dl, geometry.(geometries{k}).bt] = decsg(geometry.(geometries{k}).description, geometry.(geometries{k}).formula, geometry.(geometries{k}).names);
            
            tmp_nshapes = 1;
            tmp_description = description(:, 1);
            tmp_names = names(1, :)';
            tmp_formula = '(polygon)';
            [tmp_dl, tmp_bt] = decsg(tmp_description, tmp_formula, tmp_names);
            geometry.(geometries{k}).dl = geometry.(geometries{k}).dl;

            geometry.(geometries{k}).dl = [[tmp_dl; zeros(size(geometry.(geometries{k}).dl, 1) - 7, 4)], geometry.(geometries{k}).dl];
            geometry.(geometries{k}).dl(6 : 7, 5 : end) = geometry.(geometries{k}).dl(6 : 7, 5 : end) + 1;
        else
            geometry.(geometries{k}).n = nshapes;
            geometry.(geometries{k}).description = [description(:, 2 : end), description(:, 1)];
            geometry.(geometries{k}).names = [names(2 : end, :); names(1, :)]';
            formula = ['(', formula(14 : end), formula(1 : 10)];
            geometry.(geometries{k}).formula = append(formula, ')');
            [geometry.(geometries{k}).dl, geometry.(geometries{k}).bt] = decsg(geometry.(geometries{k}).description, geometry.(geometries{k}).formula, geometry.(geometries{k}).names);
        end
    else
        geometry.(geometries{k}).n = nshapes;
        geometry.(geometries{k}).description = description;
        geometry.(geometries{k}).names = names';
        formula = ['(', formula(4 : end)];
        geometry.(geometries{k}).formula = append(formula, ')');
        [geometry.(geometries{k}).dl, geometry.(geometries{k}).bt] = decsg(geometry.(geometries{k}).description, geometry.(geometries{k}).formula, geometry.(geometries{k}).names);
    end
    % [geometry.(geometries{k}).dl, geometry.(geometries{k}).bt] = decsg(geometry.(geometries{k}).description, geometry.(geometries{k}).formula, geometry.(geometries{k}).names);

    if 0
        model = createpde(1);
        geometryFromEdges(model, geometry.(geometries{k}).dl);
        figure
        pdegplot(geometry.(geometries{k}).dl, 'VertexLabels', 'off', 'EdgeLabels', 'off', 'FaceLabels', 'off')
        axis equal
        axis off
    end
end


% description = [];
% formula = '(';
% 
% shapes = fieldnames(geometry);
% switch parameters.lattice.lattice_type
%     case 'Square'
%         max_nrow = 10;
%     case 'Hexagonal'
%         max_nrow = 8;
% end
% max_nrow = 10;
% ngeom = 0;
% 
% for i = 1 : length(shapes)
%     shapes_i = fieldnames(geometry.(shapes{i})); 
%     for j = 1 : length(shapes_i)
%         geom = geometry.(shapes{i}).(shapes_i{j});
%         switch shapes{i}
%             case 'circle'
%                 eval_string = [shapes_i{j} ' = [1; geom.center(:); geom.radius; zeros(max_nrow - 4, 1)];'];
%             case 'polygon'
%                 eval_string = [shapes_i{j} ' = [2; size(geom.vertex, 1); geom.vertex(:); zeros(max_nrow - (2 * size(geom.vertex, 1) + 2), 1)];'];
%             case 'rectangle'
%                 eval_string = [shapes_i{j} ' = [3; 4; geom.vertex(:); zeros(max_nrow - 10, 1)];'];
%             case 'ellipse'
%                 eval_string = [shapes_i{j} ' = [4; geom.center(:); geom.semiaxes(:); geom.angle; zeros(max_nrow - 6, 1)];'];
%             case 'boundary'
%                 eval_string = [shapes_i{j} ' = [2; geom.n; geom.vertex(:); zeros(max_nrow - (2 * geom.n + 2), 1)];'];
%         end
%         eval(eval_string);
% 
%         eval(['description = [description, ' shapes_i{j} '];']);
%         if i == 1 && j == 1
%             names = char(shapes_i{j});
%         else
%             names = char(names, shapes_i{j});
%         end
%         formula = append(formula, shapes_i{j});
%         formula = append(formula, ' + ');
% 
%         ngeom = ngeom + 1;
%     end
% end
% 
% geometry.n = ngeom;
% geometry.description = description;
% geometry.names = names';
% formula(end - 2 : end) = '';
% geometry.formula = append(formula, ')');
% [geometry.dl, geometry.bt] = decsg(geometry.description, geometry.formula, geometry.names);



% %% plot geometry
% if flag.plot_geometry
%     pdegplot(geometry.dl, 'EdgeLabels', 'off', 'FaceLabels', 'off')
%     axis equal
% end






% %% define geometry
% rectangle1_designation = 3;
% rectangle1_lines_number = 4;
% rectangle1_start = [0, 0];
% rectangle1_length = [1, 1];
% circle_designation = 1;
% circle_center = [0.5, 1];
% circle_radius = 0.3;
% 
% rectangle1 = [rectangle1_designation; rectangle1_lines_number; ...
%     rectangle1_start(1); rectangle1_start(1) + rectangle1_length(1); ...
%     rectangle1_start(1) + rectangle1_length(1); rectangle1_start(1); ...
%     rectangle1_start(2); rectangle1_start(2); ...
%     rectangle1_start(2) + rectangle1_length(2); rectangle1_start(2) + rectangle1_length(2)];
% circle1 = [circle_designation; circle_center(1); circle_center(2); circle_radius];
% circle1 = [circle1; zeros(length(rectangle1) - length(circle1), 1)];
% 
% geometry.description = [rectangle1, circle1];
% geometry.names = char('rectangle1', 'circle1');
% geometry.names = geometry.names';
% geometry.formula = '(rectangle1 + circle1)';
% [geometry.dl, geometry.bt] = decsg(geometry.description, geometry.formula, geometry.names);

% 
% rectangle1 = [rectangle1_designation; rectangle1_lines_number; ...
%               rectangle1_start(1); rectangle1_start(1) + rectangle1_length(1); ...
%               rectangle1_start(1) + rectangle1_length(1); rectangle1_start(1); ...
%               rectangle1_start(2); rectangle1_start(2); ...
%               rectangle1_start(2) + rectangle1_length(2); rectangle1_start(2) + rectangle1_length(2)];
% circle1 = [circle_designation; circle_center(1); circle_center(2); circle_radius];
% circle1 = [circle1; zeros(length(rectangle1) - length(circle1), 1)];
% circle2 = [circle_designation; circle_center(1); circle_center(2); circle_radius - 0.05];
% circle2 = [circle2; zeros(length(rectangle1) - length(circle2), 1)];
% circle3 = [circle_designation; circle_center(1); circle_center(2); circle_radius + 0.05];
% circle3 = [circle3; zeros(length(rectangle1) - length(circle3), 1)];
% 
% geometry.description = [rectangle1, circle1, circle2, circle3];
% geometry.names = char('rectangle1', 'circle1', 'circle2', 'circle3');
% geometry.names = geometry.names';
% geometry.formula = '(rectangle1 + circle2) + (circle1 - circle2) + (circle3 - circle1)';
% [geometry.dl, geometry.bt] = decsg(geometry.description, geometry.formula, geometry.names);
% 
% 
% rectangle1 = [rectangle1_designation; rectangle1_lines_number; ...
%               rectangle1_start(1); rectangle1_start(1) + rectangle1_length(1); ...
%               rectangle1_start(1) + rectangle1_length(1); rectangle1_start(1); ...
%               rectangle1_start(2); rectangle1_start(2); ...
%               rectangle1_start(2) + rectangle1_length(2); rectangle1_start(2) + rectangle1_length(2)];
% circle1 = [circle_designation; circle_center(1); circle_center(2); circle_radius];
% circle1 = [circle1; zeros(length(rectangle1) - length(circle1), 1)];
% rectangle2 = rectangle1;
% rectangle2(3 : 6) = rectangle2(3 : 6) + parameter.len_a1;
% circle2 = circle1;
% circle2(2 : 3) = circle2(2 : 3) + parameter.shift;
% 
% geometry.description = [rectangle1, circle1, rectangle2, circle2];
% geometry.names = char('rectangle1', 'circle1', 'rectangle2', 'circle2');
% geometry.names = geometry.names';
% geometry.formula = '(rectangle1 + circle1) + (rectangle2 + circle2)';
% [geometry.dl, geometry.bt] = decsg(geometry.description, geometry.formula, geometry.names);

end