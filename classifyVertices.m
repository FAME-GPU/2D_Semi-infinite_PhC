function [models, geometry, meshs] = classifyVertices(models, geometry, parameters)

geometries = fieldnames(geometry);
ngeom = length(geometries);

for gg = 1 : ngeom

    geom = geometry.(geometries{gg}).model.Geometry;
    nodes = models.(geometries{gg}).model.Mesh.Nodes;
    elements = models.(geometries{gg}).model.Mesh.Elements;
    model = models.(geometries{gg}).model;
    boundary_conditions = parameters.math.(geometries{gg}).boundary_conditions;
    n_boundary = length(boundary_conditions);
    vertex = geometry.(geometries{gg}).boundary.polygon.vertex;

    % the edge must be indexed counterclockwisely
    edge1 = vertex(1 : 2, :).';
    edge2 = vertex(2 : 3, :).';
    edge3 = vertex(3 : 4, :).';
    edge4 = vertex([4 1], :).';

    %% Edges
    % indexed counterclockwisely, the first edge is 0 -> a1
    index_edge = cell(geom.NumEdges, 1);
    actual_index_edge = zeros(geom.NumEdges, 1);
    for i = 1 : geom.NumEdges
        index_edge{i} = findNodes(model.Mesh, "region", "Edge", i);
        actual_index_edge(i) = isOnBoundary(nodes(:, index_edge{i}(2)), edge1, edge2, edge3, edge4);
    end
    edge_index_for_outer_boundary = cell(n_boundary, 1);
    node_index_for_outer_boundary = cell(n_boundary, 1);
    % sort the indices of edges
    for i = 1 : n_boundary
        edge_index_for_outer_boundary{i} = find(actual_index_edge == i);
        for j = 1 : length(edge_index_for_outer_boundary{i})
            node_index_for_outer_boundary{i} = [node_index_for_outer_boundary{i}, index_edge{edge_index_for_outer_boundary{i}(j)}];
        end
    end
    % tmp_edge_index_for_outer_boundary = [];
    % for i = 1 : length(edge_index_for_outer_boundary)
    %     tmp_edge_index_for_outer_boundary = [tmp_edge_index_for_outer_boundary; edge_index_for_outer_boundary{i}];
    % end
    % edge_index_for_inner_boundary = setdiff(1 : geom.NumEdges, tmp_edge_index_for_outer_boundary);
    % node_index_for_inner_boundary = cell(geom.NumEdges - n_boundary, 1);
    % for i = 1 : geom.NumEdges - n_boundary
    %     node_index_for_inner_boundary{i} = index_edge{edge_index_for_inner_boundary(i)};
    % end 

    % quasiperiodic, find nodes near the second quasi-periodic boundary (positive direction)
    flag_quasiperiodic = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'quasiperiodic')
            flag_quasiperiodic = [flag_quasiperiodic, i];
        end
        % break; % need only the first one
    end
    if ~isempty(flag_quasiperiodic)
        flag_quasiperiodic = flip(flag_quasiperiodic);
        num_quasiperiodic = length(flag_quasiperiodic);
        node_index_quasiperiodic = cell(num_quasiperiodic, 1);
        element_index_quasiperiodic = cell(num_quasiperiodic, 1);
        node_index_near_quasiperiodic = cell(num_quasiperiodic, 1);
        for i = 1 : length(flag_quasiperiodic)
            eval(['node_index_quasiperiodic{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_quasiperiodic(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_quasiperiodic{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_quasiperiodic{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_quasiperiodic{' num2str(i) '} = unique(elements(:, element_index_quasiperiodic{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_quasiperiodic{' num2str(i) '} = unique(node_index_near_quasiperiodic{' num2str(i) '}, ''stable'').'';']);
            eval(['node_index_near_quasiperiodic{' num2str(i) '} = unique(setdiff(node_index_near_quasiperiodic{' num2str(i) '}, node_index_quasiperiodic{' num2str(i) '}, ''stable''), ''stable'');']);
        end
    end

    % Dirichlet, the same operation as that of quasiperiodic
    flag_Dirichlet = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'Dirichlet')
            flag_Dirichlet = [flag_Dirichlet; i];
        end
    end
    if ~isempty(flag_Dirichlet)
        num_Dirichlet = length(flag_Dirichlet);
        node_index_Dirichlet = cell(num_Dirichlet, 1);
        element_index_Dirichlet = cell(num_Dirichlet, 1);
        node_index_near_Dirichlet = cell(num_Dirichlet, 1);
        for i = 1 : length(flag_Dirichlet)
            eval(['node_index_Dirichlet{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_Dirichlet(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_Dirichlet{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_Dirichlet{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_Dirichlet{' num2str(i) '} = unique(elements(:, element_index_Dirichlet{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_Dirichlet{' num2str(i) '} = unique(node_index_near_Dirichlet{' num2str(i) '}, ''stable'');']);
            eval(['node_index_near_Dirichlet{' num2str(i) '} = unique(setdiff(node_index_near_Dirichlet{' num2str(i) '}, node_index_Dirichlet{' num2str(i) '}, ''stable''), ''stable'');']);
        end
    end

    % Neumann, the same operation as that of quasiperiodic
    flag_Neumann = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'Neumann')
            flag_Neumann = [flag_Neumann; i];
        end
    end
    if ~isempty(flag_Neumann)
        num_Neumann = length(flag_Neumann);
        node_index_Neumann = cell(num_Neumann, 1);
        element_index_Neumann = cell(num_Neumann, 1);
        node_index_near_Neumann = cell(num_Neumann, 1);
        for i = 1 : length(flag_Neumann)
            eval(['node_index_Neumann{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_Neumann(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_Neumann{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_Neumann{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_Neumann{' num2str(i) '} = unique(elements(:, element_index_Neumann{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_Neumann{' num2str(i) '} = unique(node_index_near_Neumann{' num2str(i) '}, ''stable'');']);
            eval(['node_index_near_Neumann{' num2str(i) '} = unique(setdiff(node_index_near_Neumann{' num2str(i) '}, node_index_Neumann{' num2str(i) '}, ''stable''), ''stable'');']);
        end
    end

    % Robin, the same operation as that of quasiperiodic
    flag_Robin = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'Robin')
            flag_Robin = [flag_Robin; i];
        end
    end
    if ~isempty(flag_Robin)
        num_Robin = length(flag_Robin);
        node_index_Robin = cell(num_Robin, 1);
        element_index_Robin = cell(num_Robin, 1);
        node_index_near_Robin = cell(num_Robin, 1);
        for i = 1 : length(flag_Robin)
            eval(['node_index_Robin{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_Robin(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_Robin{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_Robin{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_Robin{' num2str(i) '} = unique(elements(:, element_index_Robin{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_Robin{' num2str(i) '} = unique(node_index_near_Robin{' num2str(i) '}, ''stable'');']);
            eval(['node_index_near_Robin{' num2str(i) '} = unique(setdiff(node_index_near_Robin{' num2str(i) '}, node_index_Robin{' num2str(i) '}, ''stable''), ''stable'');']);
        end
    end

    % none, the same operation as that of quasiperiodic
    flag_none = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'none')
            flag_none = [flag_none; i];
        end
    end
    if ~isempty(flag_none)
        num_none = length(flag_none);
        node_index_none = cell(num_none, 1);
        element_index_none = cell(num_none, 1);
        node_index_near_none = cell(num_none, 1);
        for i = 1 : length(flag_none)
            eval(['node_index_none{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_none(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_none{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_none{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_none{' num2str(i) '} = unique(elements(:, element_index_none{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_none{' num2str(i) '} = unique(node_index_near_none{' num2str(i) '}, ''stable'');']);
            eval(['node_index_near_none{' num2str(i) '} = unique(setdiff(node_index_near_none{' num2str(i) '}, node_index_none{' num2str(i) '}, ''stable''), ''stable'');']);
        end 
    end

    % infinite, the same operation as that of quasiperiodic
    flag_infinite = [];
    for i = 1 : n_boundary
        if strcmp(boundary_conditions{i}, 'infinite')
            flag_infinite = [flag_infinite; i];
        end
    end
    if ~isempty(flag_infinite)
        num_infinite = length(flag_infinite);
        node_index_infinite = cell(num_infinite, 1);
        element_index_infinite = cell(num_infinite, 1);
        node_index_near_infinite = cell(num_infinite, 1);
        num_infinite_opposite = length(flag_infinite);
        node_index_infinite_opposite = cell(num_infinite, 1);
        element_index_infinite_opposite = cell(num_infinite, 1);
        node_index_near_infinite_opposite = cell(num_infinite, 1);
        % find the opposite boundary for the infinite boundary
        if flag_infinite > 2
            flag_infinite_opposite = flag_infinite - 2;
        else
            flag_infinite_opposite = flag_infinite + 2;
        end
        % would not be the quasi-periodic boundary condition
        opposite_boundary_condition = boundary_conditions{flag_infinite_opposite};
        if strcmp(opposite_boundary_condition, 'quasi-periodic')
            error('The opposite boundary condition for infinite can not be quasi-periodic.');
        end
        for i = 1 : length(flag_infinite)
            eval(['node_index_infinite{' num2str(i) '} = unique(node_index_for_outer_boundary{flag_infinite(' num2str(i) ')}, ''stable'');']);
            eval(['element_index_infinite{' num2str(i) '} = unique(find((sum(ismember(elements, node_index_infinite{' num2str(i) '}), 1) > 0) ~= 0), ''stable'');']);
            eval(['node_index_near_infinite{' num2str(i) '} = unique(elements(:, element_index_infinite{' num2str(i) '}), ''stable'');']);
            eval(['node_index_near_infinite{' num2str(i) '} = unique(node_index_near_infinite{' num2str(i) '}, ''stable'').'';']);
            eval(['node_index_near_infinite{' num2str(i) '} = unique(setdiff(node_index_near_infinite{' num2str(i) '}, node_index_infinite{' num2str(i) '}, ''stable''), ''stable'');']);

            eval(['node_index_infinite_opposite{' num2str(i) '} = unique(node_index_' opposite_boundary_condition '{' num2str(i) '}, ''stable'');']);
            eval(['element_index_infinite_opposite{' num2str(i) '} = unique(element_index_' opposite_boundary_condition '{' num2str(i) '}, ''stable'');']);
            eval(['node_index_near_infinite_opposite{' num2str(i) '} = unique(node_index_near_' opposite_boundary_condition '{' num2str(i) '}, ''stable'');']);
        end
    end

    % % none, the same operation as that of quasiperiodic
    % flag_none = [];
    % for i = 1 : n_boundary
    %     if strcmp(boundary_conditions{i}, 'none')
    %         flag_none = [flag_none; i];
    %     end
    % end
    % if ~isempty(flag_none)
    %     num_none = length(flag_none);
    %     node_index_none = cell(num_none, 1);
    %     element_index_none = cell(num_none, 1);
    %     node_index_near_none = cell(num_none, 1);
    %     num_none_opposite = length(flag_none);
    %     node_index_none_opposite = cell(num_none, 1);
    %     element_index_none_opposite = cell(num_none, 1);
    %     node_index_near_none_opposite = cell(num_none, 1);
    %     % find the opposite boundary for the none boundary
    %     if flag_none > 2
    %         flag_none_opposite = flag_none - 2;
    %     else
    %         flag_none_opposite = flag_none + 2;
    %     end
    %     % would not be the quasi-periodic boundary condition
    %     opposite_boundary_condition = boundary_conditions{flag_none_opposite};
    %     if strcmp(opposite_boundary_condition, 'quasi-periodic')
    %         error('The opposite boundary condition for none can not be quasi-periodic.');
    %     end
    %     for i = 1 : length(flag_none)
    %         eval(['node_index_none{' num2str(i) '} = node_index_for_outer_boundary{flag_none(' num2str(i) ')};']);
    %         eval(['element_index_none{' num2str(i) '} = find((sum(ismember(elements, node_index_none{' num2str(i) '}), 1) > 0) ~= 0);']);
    %         eval(['node_index_near_none{' num2str(i) '} = elements(:, element_index_none{' num2str(i) '});']);
    %         eval(['node_index_near_none{' num2str(i) '} = unique(node_index_near_none{' num2str(i) '}, ''stable'').'';']);
    %         eval(['node_index_near_none{' num2str(i) '} = setdiff(node_index_near_none{' num2str(i) '}, node_index_none{' num2str(i) '}, ''stable'');']);
    % 
    %         eval(['node_index_none_opposite{' num2str(i) '} = node_index_' opposite_boundary_condition '{' num2str(i) '};']);
    %         eval(['element_index_none_opposite{' num2str(i) '} = element_index_' opposite_boundary_condition '{' num2str(i) '};']);
    %         eval(['node_index_near_none_opposite{' num2str(i) '} = node_index_near_' opposite_boundary_condition '{' num2str(i) '};']);
    %     end
    % end

    

    %% Elements
    % the first index means the indices of nodes and elements located in the air
    index_1 = find(nodes(1, :) == 0);
    index_2 = find(nodes(2, :) == 0);
    index_0 = intersect(index_1, index_2);
    node_index_for_face = cell(geom.NumFaces, 1);
    flag_face = 0;
    for i = 1 : geom.NumFaces
        node_index_for_face{i} = findNodes(model.Mesh, "region", "Face", i);
        if ismember(index_0, node_index_for_face{i})
            flag_face = i;
        end
    end
    if flag_face ~= 1
        tmp = node_index_for_face{1};
        node_index_for_face{1} = node_index_for_face{flag_face};
        node_index_for_face{flag_face} = tmp;
    end

    elements = [elements; zeros(1, size(elements, 2))]; % the fourth row stores the indices of faces where the element lies in
    element_index_for_face = cell(geom.NumFaces, 1);
    flag_face = 0;
    for i = 1 : geom.NumFaces
        element_index_for_face{i} = findElements(model.Mesh, "region", "Face", i);
        elements(4, element_index_for_face{i}) = i;
        if ismember(index_0, elements(1 : 3, element_index_for_face{i}))
            flag_face = i;
        end
    end
    if flag_face ~= 1
        tmp = element_index_for_face{1};
        element_index_for_face{1} = element_index_for_face{flag_face};
        element_index_for_face{flag_face} = tmp;

        elements(4, element_index_for_face{1}) = 1;
        elements(4, element_index_for_face{flag_face}) = flag_face;
    end


    %% specially consider the current case with 'quasiperiodic', 'infinite', 'quasiperiodic', 'Dirichlet/Neumann/Robin'
    % Sort and record the required information
    [~, index_sort] = sort(nodes(1, node_index_quasiperiodic{1}), 'ascend');
    index.node_quasiperiodic1 = node_index_quasiperiodic{1}(index_sort);
    [~, index_sort] = sort(nodes(1, node_index_quasiperiodic{2}), 'ascend');
    index.node_quasiperiodic2 = node_index_quasiperiodic{2}(index_sort);
    [~, index_sort] = sort(nodes(1, node_index_near_quasiperiodic{1}), 'ascend');
    index.node_near_quasiperiodic1 = node_index_near_quasiperiodic{1}(index_sort);
    [~, index_sort] = sort(nodes(1, node_index_near_quasiperiodic{2}), 'ascend');
    index.node_near_quasiperiodic2 = node_index_near_quasiperiodic{2}(index_sort);
    index.node_near_quasiperiodic2_plus_quasiperiodic1 = [index.node_near_quasiperiodic2, index.node_quasiperiodic1];
    index.node_near_quasiperiodic2_plus_quasiperiodic2 = [index.node_near_quasiperiodic2, index.node_quasiperiodic2];

    if ~isempty(flag_infinite)
        [~, index_sort] = sort(nodes(1, node_index_infinite{1}), 'ascend');
        index.node_infinite = node_index_infinite{1}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_near_infinite{1}), 'ascend');
        index.node_near_infinite = node_index_near_infinite{1}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_infinite_opposite{1}), 'ascend');
        index.node_infinite_opposite = node_index_infinite_opposite{1}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_near_infinite_opposite{1}), 'ascend');
        index.node_near_infinite_opposite = node_index_near_infinite_opposite{1}(index_sort);
        index.element_infinite_opposite = element_index_infinite_opposite{1};
    
        index.node_purely_infinite = setdiff(index.node_infinite, index.node_quasiperiodic2, 'stable');
        index.node_purely_near_infinite = setdiff(index.node_near_infinite, index.node_quasiperiodic2, 'stable');
        index.node_purely_infinite_opposite = setdiff(index.node_infinite_opposite, index.node_quasiperiodic2, 'stable');
        index.node_purely_near_infinite_opposite = setdiff(index.node_near_infinite_opposite, index.node_quasiperiodic2, 'stable');
        if size(index.node_purely_near_infinite_opposite, 1) > 1
            index.node_purely_near_infinite_opposite = index.node_purely_near_infinite_opposite';
        end
        index.node_purely_near_infinite_plus_infinite = [index.node_purely_near_infinite, index.node_purely_infinite];
    
        index.free = 1 : size(nodes, 2);
        index.restricted = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                            index.node_quasiperiodic2, index.node_purely_near_infinite, index.node_purely_infinite];
        index.free = setdiff(index.free, index.restricted);
    end
    if length(flag_none) == 2
        [~, index_sort] = sort(nodes(1, node_index_none{1}), 'ascend');
        index.node_none{1} = node_index_none{1}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_near_none{1}), 'ascend');
        index.node_near_none{1} = node_index_near_none{1}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_none{2}), 'ascend');
        index.node_none{2} = node_index_none{2}(index_sort);
        [~, index_sort] = sort(nodes(1, node_index_near_none{2}), 'ascend');
        index.node_near_none{2} = node_index_near_none{2}(index_sort);
        index.element_none{1} = element_index_none{1};
        index.element_none{2} = element_index_none{2};
    
        index.node_purely_none{1} = setdiff(index.node_none{1}, index.node_quasiperiodic2, 'stable');
        index.node_purely_near_none{1} = setdiff(index.node_near_none{1}, index.node_quasiperiodic2, 'stable');
        index.node_purely_none{2} = setdiff(index.node_none{2}, index.node_quasiperiodic2, 'stable');
        index.node_purely_near_none{2} = setdiff(index.node_near_none{2}, index.node_quasiperiodic2, 'stable');
        if size(index.node_purely_near_none{1}, 1) > 1
            index.node_purely_near_none{1} = index.node_purely_near_none{1}';
        end
        if size(index.node_purely_near_none{2}, 1) > 1
            index.node_purely_near_none{2} = index.node_purely_near_none{2}';
        end
        if size(index.node_near_none{1}, 1) > 1
            index.node_near_none{1} = index.node_near_none{1}';
        end
        index.node_purely_near_none_plus_none{1} = [index.node_purely_near_none{1}, index.node_purely_none{1}];
        index.node_purely_near_none_plus_none{2} = [index.node_purely_near_none{2}, index.node_purely_none{2}];
    
        index.free = 1 : size(nodes, 2);
        index.restricted = [index.node_purely_none{2}, index.node_purely_near_none{2}, ... 
                            index.node_quasiperiodic2, index.node_purely_near_none{1}, index.node_purely_none{1}];
        index.free = setdiff(index.free, index.restricted);
    end

    if ngeom == 1
        if ~isempty(flag_Dirichlet)
            index.free_end = [index.node_purely_near_infinite_opposite, index.free, index.node_purely_near_infinite];
            index.free_begin = [index.free_end, index.node_purely_infinite];
            index.free = index.free_begin;
            % index.free_begin = [index.node_purely_near_infinite_opposite, index.free, index.node_purely_near_infinite];
            % index.free = [index.node_purely_infinite_opposite, index.free_begin];
            % index.free_end = index.free;
        elseif ~isempty(flag_Neumann)
            index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, index.free, index.node_purely_near_infinite];
            index.free_begin = index.free;
            index.free_end = [index.free, index.node_purely_infinite];
        else
            index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, index.free, index.node_purely_near_infinite];
            index.free_begin = index.free;
            index.free_end = [index.free, index.node_purely_infinite];
        end
    elseif ngeom == 2
        if strcmp(parameters.math.equation, 'TM')
            flag_Dirichlet = flag_infinite;
        elseif strcmp(parameters.math.equation, 'TE')
            flag_Neumann = flag_infinite;
        else

        end
        if gg == 1
            if ~isempty(flag_Dirichlet)
                index.free_end = [index.node_purely_near_infinite, index.free, index.node_purely_near_infinite_opposite];
                index.free = [index.node_purely_infinite, index.free_end];
                index.free_begin = index.free;
            elseif ~isempty(flag_Neumann)
                index.free_end = [index.node_purely_infinite, index.node_purely_near_infinite, index.free, ... 
                                  index.node_purely_near_infinite_opposite];
                index.free = index.free_end;
                index.free_begin = index.free;
            else
                index.free_end = [index.node_purely_infinite, index.node_purely_near_infinite, index.free, ... 
                                  index.node_purely_near_infinite_opposite];
                index.free = index.free_end;
                index.free_begin = index.free;
            end
        elseif gg == 2
            if ~isempty(flag_Dirichlet)
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = index.free;
            elseif ~isempty(flag_Neumann)
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = [index.free, index.node_purely_infinite];
            else
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = [index.free, index.node_purely_infinite];
            end
        end
    elseif ngeom == 3
        if strcmp(parameters.math.equation, 'TM')
            flag_Dirichlet = flag_infinite;
        elseif strcmp(parameters.math.equation, 'TE')
            flag_Neumann = flag_infinite;
        else

        end
        if gg == 1
            if ~isempty(flag_Dirichlet)
                index.free = [index.node_purely_near_infinite, index.free, index.node_purely_near_infinite_opposite, ... 
                              index.node_purely_infinite_opposite];
                index.free_begin = index.free;
                index.free_end = index.free;
            elseif ~isempty(flag_Neumann)
                index.free_begin = [index.node_purely_near_infinite, index.free, ... 
                                  index.node_purely_near_infinite_opposite, index.node_purely_infinite_opposite];
                index.free = index.free_begin;
                index.free_end = [index.node_purely_infinite, index.free];
            else
                index.free_begin = [index.node_purely_near_infinite, index.free, ... 
                                  index.node_purely_near_infinite_opposite, index.node_purely_infinite_opposite];
                index.free = index.free_begin;
                index.free_end = [index.node_purely_infinite, index.free];
            end
        elseif gg == 2  % K_C use all DOF in geometry 2 except for nodes of boundary condition 'none'
            index.free = [index.node_purely_near_none{1}, index.free, index.node_purely_near_none{2}];
            index.free_begin = index.free;
            index.free_end = index.free;
        elseif gg == 3
            if ~isempty(flag_Dirichlet)
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = index.free;
            elseif ~isempty(flag_Neumann)
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = [index.free, index.node_purely_infinite];
            else
                index.free = [index.node_purely_infinite_opposite, index.node_purely_near_infinite_opposite, ... 
                              index.free, index.node_purely_near_infinite];
                index.free_begin = index.free;
                index.free_end = [index.free, index.node_purely_infinite];
            end
        end
    end

    tmp_ismember = ismember(elements(1 : 3, element_index_quasiperiodic{1}), index.node_quasiperiodic1);
    index.quasiperiodic1_nodes_in_elements = zeros(3, size(elements, 2));
    index.quasiperiodic1_nodes_in_elements(:, element_index_quasiperiodic{1}) = tmp_ismember;
    index.quasiperiodic1_nodes_in_elements = logical(index.quasiperiodic1_nodes_in_elements);
    tmp_ismember = ismember(elements(1 : 3, element_index_quasiperiodic{2}), index.node_quasiperiodic2);
    index.quasiperiodic2_nodes_in_elements = zeros(3, size(elements, 2));
    index.quasiperiodic2_nodes_in_elements(:, element_index_quasiperiodic{2}) = tmp_ismember;
    index.quasiperiodic2_nodes_in_elements = logical(index.quasiperiodic2_nodes_in_elements);
    tmp_ismember = ismember(elements(1 : 3, element_index_infinite{1}), index.node_infinite);
    index.infinite_nodes_in_elements = zeros(3, size(elements, 2));
    index.infinite_nodes_in_elements(:, element_index_infinite{1}) = tmp_ismember;
    index.infinite_nodes_in_elements = logical(index.infinite_nodes_in_elements);
    tmp_ismember = ismember(elements(1 : 3, element_index_infinite_opposite{1}), index.node_infinite_opposite);
    index.infinite_opposite_nodes_in_elements = zeros(3, size(elements, 2));
    index.infinite_opposite_nodes_in_elements(:, element_index_infinite_opposite{1}) = tmp_ismember;
    index.infinite_opposite_nodes_in_elements = logical(index.infinite_opposite_nodes_in_elements);

    % find(index.quasiperiodic2_nodes_in_elements + index.infinite_opposite_nodes_in_elements == 2)


    mesh.nodes = nodes;
    mesh.elements = elements;
    mesh.elements_simple = elements(1 : 3, :);
    mesh.elements_infinite_opposite = elements(:, index.element_infinite_opposite);
    mesh.elements_infinite_opposite_simple = mesh.elements_infinite_opposite(1 : 3, :);
    mesh.index = index;
    mesh.MaxElementSize = model.Mesh.MaxElementSize;
    mesh.MinElementSize = model.Mesh.MinElementSize;
    mesh.GeometricOrder = model.Mesh.GeometricOrder;

    mesh.free_nodes = nodes(:, index.free);
    mesh.free_begin_nodes = nodes(:, index.free_begin);
    mesh.free_end_nodes = nodes(:, index.free_end);

    mesh.n_nodes = size(nodes, 2);
    mesh.n_elements = size(elements, 2);
    mesh.n_infinite_opposite_elements = size(mesh.elements_infinite_opposite, 2);
    mesh.n_free_nodes = size(mesh.free_nodes, 2);
    mesh.n_free_begin_nodes = size(mesh.free_begin_nodes, 2);
    mesh.n_free_end_nodes = size(mesh.free_end_nodes, 2);

    % mesh.n_free_elements = size(mesh.free_elements, 2);
    % mesh.n_purely_interior_elements = length(index.element_purely_interior);
    % mesh.n_purely_infinite_elements = length(index.element_purely_infinite);
    % mesh.n_purely_infinite_opposite_elements = length(index.element_purely_infinite_opposite);
    % mesh.n_quasiperiodic_elements = zeros(2, 1);
    % mesh.n_quasiperiodic_elements(1) = length(index.element_quasiperiodic{1});
    % mesh.n_quasiperiodic_elements(2) = length(index.element_quasiperiodic{2});

    % ppp = 1;


    meshs.(geometries{gg}) = mesh;
    meshs.(geometries{gg}).index = index;

end




















% geom = model.Geometry;
% nodes = model.Mesh.Nodes;
% elements = model.Mesh.Elements;
% boundary_conditions = parameters.math.boundary_condition.self;
% n_boundary = length(boundary_conditions);
% vertex = parameters.geometry.boundary.polygon.vertex;
% 
% % the edge must be indexed counterclockwisely
% edge1 = vertex(1 : 2, :).';
% edge2 = vertex(2 : 3, :).';
% edge3 = vertex(3 : 4, :).';
% edge4 = vertex([4 1], :).';
% 
% %% Edges
% % indexed counterclockwisely, the first edge is 0 -> a1
% index_edge = cell(geom.NumEdges, 1);
% actual_index_edge = zeros(geom.NumEdges, 1);
% for i = 1 : geom.NumEdges
%     index_edge{i} = findNodes(model.Mesh, "region", "Edge", i);
%     actual_index_edge(i) = isOnBoundary(nodes(:, index_edge{i}(2)), edge1, edge2, edge3, edge4);
% end
% edge_index_for_outer_boundary = zeros(n_boundary, 1);
% node_index_for_outer_boundary = cell(n_boundary, 1); 
% % sort the indices of edges
% for i = 1 : n_boundary
%     edge_index_for_outer_boundary(i) = find(actual_index_edge == i);
%     node_index_for_outer_boundary{i} = index_edge{edge_index_for_outer_boundary(i)};
% end
% edge_index_for_inner_boundary = setdiff(1 : geom.NumEdges, edge_index_for_outer_boundary);
% node_index_for_inner_boundary = cell(geom.NumEdges - n_boundary, 1);
% for i = 1 : geom.NumEdges - n_boundary
%     node_index_for_inner_boundary{i} = index_edge{edge_index_for_inner_boundary(i)};
% end
% 
% % quasiperiodic, find nodes near the second quasi-periodic boundary (positive direction)
% flag_quasiperiodic = [];
% for i = 1 : n_boundary
%     if strcmp(boundary_conditions{i}, 'quasiperiodic')
%         flag_quasiperiodic = [flag_quasiperiodic, i]; 
%     end
%     % break; % need only the first one
% end
% if ~isempty(flag_quasiperiodic)
%     num_quasiperiodic = length(flag_quasiperiodic);
%     node_index_quasiperiodic = cell(num_quasiperiodic, 1);
%     element_index_quasiperiodic = cell(num_quasiperiodic, 1);
%     node_index_near_quasiperiodic = cell(num_quasiperiodic, 1);
%     for i = 1 : length(flag_quasiperiodic)
%         eval(['node_index_quasiperiodic{' num2str(i) '} = node_index_for_outer_boundary{flag_quasiperiodic(' num2str(i) ')};']);
%         eval(['element_index_quasiperiodic{' num2str(i) '} = find((sum(ismember(elements, node_index_quasiperiodic{' num2str(i) '}), 1) > 0) ~= 0);']);
%         eval(['node_index_near_quasiperiodic{' num2str(i) '} = elements(:, element_index_quasiperiodic{' num2str(i) '});']);
%         eval(['node_index_near_quasiperiodic{' num2str(i) '} = unique(node_index_near_quasiperiodic{' num2str(i) '}, ''stable'').'';']);
%         eval(['node_index_near_quasiperiodic{' num2str(i) '} = setdiff(node_index_near_quasiperiodic{' num2str(i) '}, node_index_quasiperiodic{' num2str(i) '}, ''stable'');']);
%     end
% end
% 
% % Dirichlet, the same operation as that of quasiperiodic
% flag_Dirichlet = [];
% for i = 1 : n_boundary
%     if strcmp(boundary_conditions{i}, 'Dirichlet')
%         flag_Dirichlet = [flag_Dirichlet; i];
%     end
% end
% if ~isempty(flag_Dirichlet)
%     num_Dirichlet = length(flag_Dirichlet);
%     node_index_Dirichlet = cell(num_Dirichlet, 1);
%     element_index_Dirichlet = cell(num_Dirichlet, 1);
%     node_index_near_Dirichlet = cell(num_Dirichlet, 1);
%     for i = 1 : length(flag_Dirichlet)
%         eval(['node_index_Dirichlet{' num2str(i) '} = node_index_for_outer_boundary{flag_Dirichlet(' num2str(i) ')};']);
%         eval(['element_index_Dirichlet{' num2str(i) '} = find((sum(ismember(elements, node_index_Dirichlet{' num2str(i) '}), 1) > 0) ~= 0);']);
%         eval(['node_index_near_Dirichlet{' num2str(i) '} = elements(:, element_index_Dirichlet{' num2str(i) '});']);
%         eval(['node_index_near_Dirichlet{' num2str(i) '} = unique(node_index_near_Dirichlet{' num2str(i) '}, ''stable'');']);
%         eval(['node_index_near_Dirichlet{' num2str(i) '} = setdiff(node_index_near_Dirichlet{' num2str(i) '}, node_index_Dirichlet{' num2str(i) '}, ''stable'');']);
%     end
% end
% 
% % Neumann, the same operation as that of quasiperiodic
% flag_Neumann = [];
% for i = 1 : n_boundary
%     if strcmp(boundary_conditions{i}, 'Neumann')
%         flag_Neumann = [flag_Neumann; i];
%     end
% end
% if ~isempty(flag_Neumann)
%     num_Neumann = length(flag_Neumann);
%     node_index_Neumann = cell(num_Neumann, 1);
%     element_index_Neumann = cell(num_Neumann, 1);
%     node_index_near_Neumann = cell(num_Neumann, 1);
%     for i = 1 : length(flag_Neumann)
%         eval(['node_index_Neumann{' num2str(i) '} = node_index_for_outer_boundary{flag_Neumann(' num2str(i) ')};']);
%         eval(['element_index_Neumann{' num2str(i) '} = find((sum(ismember(elements, node_index_Neumann{' num2str(i) '}), 1) > 0) ~= 0);']);
%         eval(['node_index_near_Neumann{' num2str(i) '} = elements(:, element_index_Neumann{' num2str(i) '});']);
%         eval(['node_index_near_Neumann{' num2str(i) '} = unique(node_index_near_Neumann{' num2str(i) '}, ''stable'');']);
%         eval(['node_index_near_Neumann{' num2str(i) '} = setdiff(node_index_near_Neumann{' num2str(i) '}, node_index_Neumann{' num2str(i) '}, ''stable'');']);
%     end
% end
% 
% % Robin, the same operation as that of quasiperiodic
% flag_Robin = [];
% for i = 1 : n_boundary
%     if strcmp(boundary_conditions{i}, 'Robin')
%         flag_Robin = [flag_Robin; i];
%     end
% end
% if ~isempty(flag_Robin)
%     num_Robin = length(flag_Robin);
%     node_index_Robin = cell(num_Robin, 1);
%     element_index_Robin = cell(num_Robin, 1);
%     node_index_near_Robin = cell(num_Robin, 1);
%     for i = 1 : length(flag_Robin)
%         eval(['node_index_Robin{' num2str(i) '} = node_index_for_outer_boundary{flag_Robin(' num2str(i) ')};']);
%         eval(['element_index_Robin{' num2str(i) '} = find((sum(ismember(elements, node_index_Robin{' num2str(i) '}), 1) > 0) ~= 0);']);
%         eval(['node_index_near_Robin{' num2str(i) '} = elements(:, element_index_Robin{' num2str(i) '});']);
%         eval(['node_index_near_Robin{' num2str(i) '} = unique(node_index_near_Robin{' num2str(i) '}, ''stable'');']);
%         eval(['node_index_near_Robin{' num2str(i) '} = setdiff(node_index_near_Robin{' num2str(i) '}, node_index_Robin{' num2str(i) '}, ''stable'');']);
%     end
% end
% 
% % infinite, the same operation as that of quasiperiodic
% flag_infinite = [];
% for i = 1 : n_boundary
%     if strcmp(boundary_conditions{i}, 'infinite')
%         flag_infinite = [flag_infinite; i];
%     end
% end
% if ~isempty(flag_infinite)
%     num_infinite = length(flag_infinite);
%     node_index_infinite = cell(num_infinite, 1);
%     element_index_infinite = cell(num_infinite, 1);
%     node_index_near_infinite = cell(num_infinite, 1);
%     num_infinite_opposite = length(flag_infinite);
%     node_index_infinite_opposite = cell(num_infinite, 1);
%     element_index_infinite_opposite = cell(num_infinite, 1);
%     node_index_near_infinite_opposite = cell(num_infinite, 1);
%     % find the opposite boundary for the infinite boundary
%     if flag_infinite > 2
%         flag_infinite_opposite = flag_infinite - 2;
%     else
%         flag_infinite_opposite = flag_infinite + 2;
%     end
%     % would not be the quasi-periodic boundary condition
%     opposite_boundary_condition = boundary_conditions{flag_infinite_opposite};
%     if strcmp(opposite_boundary_condition, 'quasi-periodic')
%         error('The opposite boundary consition for infinite can not be quasi-periodic.');
%     end
%     for i = 1 : length(flag_infinite)
%         eval(['node_index_infinite{' num2str(i) '} = node_index_for_outer_boundary{flag_infinite(' num2str(i) ')};']);
%         eval(['element_index_infinite{' num2str(i) '} = find((sum(ismember(elements, node_index_infinite{' num2str(i) '}), 1) > 0) ~= 0);']);
%         eval(['node_index_near_infinite{' num2str(i) '} = elements(:, element_index_infinite{' num2str(i) '});']);
%         eval(['node_index_near_infinite{' num2str(i) '} = unique(node_index_near_infinite{' num2str(i) '}, ''stable'').'';']);
%         eval(['node_index_near_infinite{' num2str(i) '} = setdiff(node_index_near_infinite{' num2str(i) '}, node_index_infinite{' num2str(i) '}, ''stable'');']);
% 
%         eval(['node_index_infinite_opposite{' num2str(i) '} = node_index_' opposite_boundary_condition '{' num2str(i) '};']);
%         eval(['element_index_infinite_opposite{' num2str(i) '} = element_index_' opposite_boundary_condition '{' num2str(i) '};']);
%         eval(['node_index_near_infinite_opposite{' num2str(i) '} = node_index_near_' opposite_boundary_condition '{' num2str(i) '};']);
%     end
% end
% 
% 
% %% Elements
% % the first index means the indices of nodes and elements located in the air
% index_1 = find(nodes(1, :) == 0);
% index_2 = find(nodes(2, :) == 0);
% index_0 = intersect(index_1, index_2);
% node_index_for_face = cell(geom.NumFaces, 1);
% flag_face = 0;
% for i = 1 : geom.NumFaces
%     node_index_for_face{i} = findNodes(model.Mesh, "region", "Face", i);
%     if ismember(index_0, node_index_for_face{i})
%         flag_face = i;
%     end
% end
% if flag_face ~= 1
%     tmp = node_index_for_face{1};
%     node_index_for_face{1} = node_index_for_face{flag_face};
%     node_index_for_face{flag_face} = tmp;
% end
% 
% elements = [elements; zeros(1, size(elements, 2))]; % the fourth row stores the indices of faces where the element lies in
% element_index_for_face = cell(geom.NumFaces, 1);
% flag_face = 0;
% for i = 1 : geom.NumFaces
%     element_index_for_face{i} = findElements(model.Mesh, "region", "Face", i);
%     elements(4, element_index_for_face{i}) = i;
%     if ismember(index_0, elements(1 : 3, element_index_for_face{i}))
%         flag_face = i;
%     end
% end
% if flag_face ~= 1
%     tmp = element_index_for_face{1};
%     element_index_for_face{1} = node_index_for_face{flag_face};
%     element_index_for_face{flag_face} = tmp;
% 
%     elements(4, element_index_for_face{1}) = flag_face;
%     elements(4, element_index_for_face{flag_face}) = 1;
% end
% 
% 
% %% specially consider the current case with 'quasiperiodic', 'infinite', 'quasiperiodic', 'Dirichlet/Neumann/Robin'
% % Sort and record the required information
% [~, index_sort] = sort(nodes(1, node_index_quasiperiodic{1}), 'ascend');
% index.node_quasiperiodic1 = node_index_quasiperiodic{1}(index_sort);
% [~, index_sort] = sort(nodes(1, node_index_quasiperiodic{2}), 'ascend');
% index.node_quasiperiodic2 = node_index_quasiperiodic{2}(index_sort);
% [~, index_sort] = sort(nodes(1, node_index_near_quasiperiodic{1}), 'ascend');
% index.node_near_quasiperiodic1 = node_index_near_quasiperiodic{1}(index_sort);
% [~, index_sort] = sort(nodes(1, node_index_near_quasiperiodic{2}), 'ascend');
% index.node_near_quasiperiodic2 = node_index_near_quasiperiodic{2}(index_sort);
% index.node_near_quasiperiodic2_plus_quasiperiodic1 = [index.node_near_quasiperiodic2, index.node_quasiperiodic1];
% index.node_near_quasiperiodic2_plus_quasiperiodic2 = [index.node_near_quasiperiodic2, index.node_quasiperiodic2];
% 
% [~, index_sort] = sort(nodes(2, node_index_infinite{1}), 'ascend');
% index.node_infinite = node_index_infinite{1}(index_sort);
% [~, index_sort] = sort(nodes(2, node_index_near_infinite{1}), 'ascend');
% index.node_near_infinite = node_index_near_infinite{1}(index_sort);
% [~, index_sort] = sort(nodes(2, node_index_infinite_opposite{1}), 'ascend');
% index.node_infinite_opposite = node_index_infinite_opposite{1}(index_sort);
% [~, index_sort] = sort(nodes(2, node_index_near_infinite_opposite{1}), 'ascend');
% index.node_near_infinite_opposite = node_index_near_infinite_opposite{1}(index_sort);
% index.element_infinite_opposite = element_index_infinite_opposite{1};
% 
% index.node_purely_infinite = setdiff(index.node_infinite, index.node_quasiperiodic2, 'stable');
% index.node_purely_near_infinite = setdiff(index.node_near_infinite, index.node_quasiperiodic2, 'stable');
% index.node_purely_infinite_opposite = setdiff(index.node_infinite_opposite, index.node_quasiperiodic2, 'stable');
% index.node_purely_near_infinite_opposite = setdiff(index.node_near_infinite_opposite, index.node_quasiperiodic2, 'stable');
% index.node_purely_near_infinite_plus_infinite = [index.node_purely_near_infinite, index.node_purely_infinite];
% 
% index.free = 1 : size(nodes, 2);
% index.restricted = [index.node_purely_infinite_opposite, index.node_quasiperiodic2, index.node_near_infinite, index.node_purely_infinite];
% index.free = setdiff(index.free, index.restricted);
% if ismember('Dirichlet', parameters.math.boundary_condition.self)
%     index.free_0 = [index.free, index.node_purely_near_infinite];
%     index.free = [index.node_purely_infinite_opposite, index.free_0];
%     index.free_end = [index.free];  % , index.node_purely_infinite
% else
%     index.free_0 = [index.node_purely_infinite_opposite, index.free, index.node_purely_near_infinite];
%     index.free = index.free_0;
%     index.free_end = [index.free, index.node_purely_infinite];
% end
% 
% tmp_ismember = ismember(elements(1 : 3, element_index_quasiperiodic{1}), index.node_quasiperiodic1);
% index.quasiperiodic1_nodes_in_elements = zeros(3, size(elements, 2));
% index.quasiperiodic1_nodes_in_elements(:, element_index_quasiperiodic{1}) = tmp_ismember;
% index.quasiperiodic1_nodes_in_elements = logical(index.quasiperiodic1_nodes_in_elements);
% tmp_ismember = ismember(elements(1 : 3, element_index_quasiperiodic{2}), index.node_quasiperiodic2);
% index.quasiperiodic2_nodes_in_elements = zeros(3, size(elements, 2));
% index.quasiperiodic2_nodes_in_elements(:, element_index_quasiperiodic{2}) = tmp_ismember;
% index.quasiperiodic2_nodes_in_elements = logical(index.quasiperiodic2_nodes_in_elements);
% tmp_ismember = ismember(elements(1 : 3, element_index_infinite{1}), index.node_infinite);
% index.infinite_nodes_in_elements = zeros(3, size(elements, 2));
% index.infinite_nodes_in_elements(:, element_index_infinite{1}) = tmp_ismember;
% index.infinite_nodes_in_elements = logical(index.infinite_nodes_in_elements);
% tmp_ismember = ismember(elements(1 : 3, element_index_infinite_opposite{1}), index.node_infinite_opposite);
% index.infinite_opposite_nodes_in_elements = zeros(3, size(elements, 2));
% index.infinite_opposite_nodes_in_elements(:, element_index_infinite_opposite{1}) = tmp_ismember;
% index.infinite_opposite_nodes_in_elements = logical(index.infinite_opposite_nodes_in_elements);
% 
% % find(index.quasiperiodic2_nodes_in_elements + index.infinite_opposite_nodes_in_elements == 2)
% 
% 
% mesh.nodes = nodes;
% mesh.elements = elements;
% mesh.elements_simple = elements(1 : 3, :);
% mesh.elements_infinite_opposite = elements(:, index.element_infinite_opposite);
% mesh.elements_infinite_opposite_simple = mesh.elements_infinite_opposite(1 : 3, :);
% mesh.index = index;
% mesh.MaxElementSize = model.Mesh.MaxElementSize;
% mesh.MinElementSize = model.Mesh.MinElementSize;
% mesh.GeometricOrder = model.Mesh.GeometricOrder;
% 
% mesh.free0_nodes = nodes(:, index.free_0);
% mesh.free_nodes = nodes(:, index.free);
% 
% mesh.n_nodes = size(nodes, 2);
% mesh.n_elements = size(elements, 2);
% mesh.n_infinite_opposite_elements = size(mesh.elements_infinite_opposite, 2);
% mesh.n_free_nodes = size(mesh.free_nodes, 2);
% mesh.n_free0_nodes = size(mesh.free0_nodes, 2);
% % mesh.n_free_elements = size(mesh.free_elements, 2);
% % mesh.n_purely_interior_elements = length(index.element_purely_interior);
% % mesh.n_purely_infinite_elements = length(index.element_purely_infinite);
% % mesh.n_purely_infinite_opposite_elements = length(index.element_purely_infinite_opposite);
% % mesh.n_quasiperiodic_elements = zeros(2, 1);
% % mesh.n_quasiperiodic_elements(1) = length(index.element_quasiperiodic{1});
% % mesh.n_quasiperiodic_elements(2) = length(index.element_quasiperiodic{2});
% 
% % ppp = 1;


% index.node_for_inner_boundary = node_index_for_inner_boundary;
% index.node_for_outer_boundary = node_index_for_outer_boundary;
% if ~isempty(flag_quasiperiodic)
%     index.node_quasiperiodic = node_index_quasiperiodic;
%     index.node_near_quasiperiodic = node_index_near_quasiperiodic;
%     index.node_quasiperiodic_plus_near_quasiperiodic = cell(length(flag_quasiperiodic));
%     for i = 1 : length(flag_quasiperiodic)
%         for j = 1 : length(flag_quasiperiodic)
%             index.node_quasiperiodic_plus_near_quasiperiodic{i, j} = union(index.node_quasiperiodic{i}, ... 
%                                                                            index.node_near_quasiperiodic{j});
%         end
%     end
%     index.element_quasiperiodic = element_index_quasiperiodic;
% end
% if ~isempty(flag_Dirichlet)
%     index.node_Dirichlet = node_index_Dirichlet;
%     index.node_near_Dirichlet = node_index_near_Dirichlet;
%     index.element_Dirichlet = element_index_Dirichlet;
% end
% if ~isempty(flag_Neumann)
%     index.node_Neumann = node_index_Neumann;
%     index.node_near_Neumann = node_index_near_Neumann;
%     index.element_Neumann = element_index_Neumann;
% end
% if ~isempty(flag_Robin)
%     index.node_Robin = node_index_Robin;
%     index.node_near_Robin = node_index_near_Robin;
%     index.element_Robin = element_index_Robin;
% end
% if ~isempty(flag_infinite)
%     index.node_infinite = node_index_infinite;
%     index.node_near_infinite = node_index_near_infinite;
%     index.node_near_infinite_plus_infinite = union(index.node_infinite{1}, index.node_near_infinite{1});
%     index.element_infinite = element_index_infinite;
%     index.node_infinite_opposite = node_index_infinite_opposite;
%     index.node_near_infinite_opposite = node_index_near_infinite_opposite;
%     index.element_infinite_opposite = element_index_infinite_opposite;
% end
% index.node_for_face = node_index_for_face;
% index.element_for_face = element_index_for_face;
% 
% %% specially consider the current case with 'quasiperiodic', 'infinite', 'quasiperiodic', 'Dirichlet/Neumann/Robin'
% index.node_all_boundary = union(index.node_quasiperiodic{1}, index.node_infinite{1});
% index.node_all_boundary = union(index.node_all_boundary, index.node_quasiperiodic{2});
% index.element_all_boundary = union(index.element_quasiperiodic{1}, index.element_infinite{1});
% index.element_all_boundary = union(index.element_all_boundary, index.element_quasiperiodic{2});
% index.node_corner = cell(n_boundary, 1);
% index.node_corner{1} = intersect(index.node_quasiperiodic{1}, index.node_infinite{1});
% index.node_corner{2} = intersect(index.node_quasiperiodic{2}, index.node_infinite{1});
% index.element_corner = cell(n_boundary, 1);
% index.element_corner{1} = intersect(index.element_quasiperiodic{1}, index.element_infinite{1});
% index.element_corner{2} = intersect(index.element_quasiperiodic{2}, index.element_infinite{1});
% if ~isempty(flag_Dirichlet)
%     index.node_all_boundary = union(index.node_all_boundary, index.node_Dirichlet{1});
%     index.element_all_boundary = union(index.element_all_boundary, index.element_Dirichlet{1});
%     index.node_corner{3} = intersect(index.node_quasiperiodic{2}, index.node_Dirichlet{1});
%     index.node_corner{4} = intersect(index.node_quasiperiodic{1}, index.node_Dirichlet{1});
%     index.node_purely_infinite_opposite = setdiff(index.node_Dirichlet{1}, ... 
%                                             union(index.node_corner{3}, index.node_corner{4}));
%     index.element_corner{3} = intersect(index.element_quasiperiodic{2}, index.element_Dirichlet{1});
%     index.element_corner{4} = intersect(index.element_quasiperiodic{1}, index.element_Dirichlet{1});
%     index.element_purely_infinite_opposite = setdiff(index.element_Dirichlet{1}, ... 
%                                             union(index.element_corner{3}, index.element_corner{4}));
%     index.infinite_opposite_nodes_in_elements = ismember(elements(1 : 3, :), index.node_Dirichlet{1});
% elseif ~isempty(flag_Neumann)
%     index.node_all_boundary = union(index.node_all_boundary, index.node_Neumann{1});
%     index.element_all_boundary = union(index.element_all_boundary, index.element_Neumann{1});
%     index.node_corner{3} = intersect(index.node_quasiperiodic{2}, index.node_Neumann{1});
%     index.node_corner{4} = intersect(index.node_quasiperiodic{1}, index.node_Neumann{1});
%     index.node_purely_infinite_opposite = setdiff(index.node_Neumann{1}, ... 
%                                             union(index.node_corner{3}, index.node_corner{4}));
%     index.element_corner{3} = intersect(index.element_quasiperiodic{2}, index.element_Neumann{1});
%     index.element_corner{4} = intersect(index.element_quasiperiodic{1}, index.element_Neumann{1});
%     index.element_purely_infinite_opposite = setdiff(index.element_Neumann{1}, ... 
%                                             union(index.element_corner{3}, index.element_corner{4}));
%     index.infinite_opposite_nodes_in_elements = ismember(elements(1 : 3, :), index.node_Neumann{1});
% elseif ~isempty(flag_Robin)
%     index.node_all_boundary = union(index.node_all_boundary, index.node_Robin{1});
%     index.element_all_boundary = union(index.element_all_boundary, index.element_Robin{1});
%     index.node_corner{3} = intersect(index.node_quasiperiodic{2}, index.node_Robin{1});
%     index.node_corner{4} = intersect(index.node_quasiperiodic{1}, index.node_Robin{1});
%     index.node_purely_infinite_opposite = setdiff(index.node_Robin{1}, ... 
%                                             union(index.node_corner{3}, index.node_corner{4}));
%     index.element_corner{3} = intersect(index.element_quasiperiodic{2}, index.element_Robin{1});
%     index.element_corner{4} = intersect(index.element_quasiperiodic{1}, index.element_Robin{1});
%     index.element_purely_infinite_opposite = setdiff(index.element_Robin{1}, ... 
%                                             union(index.element_corner{3}, index.element_corner{4}));
%     index.infinite_opposite_nodes_in_elements = ismember(elements(1 : 3, :), index.node_Robin{1});
% end
% index.infinite_nodes_in_elements = ismember(elements(1 : 3, :), index.node_infinite{1});
% index.quasiperiodic_nodes_in_elements = ismember(elements(1 : 3, :), index.node_quasiperiodic{1});
% index.quasiperiodic2_nodes_in_elements = ismember(elements(1 : 3, :), index.node_quasiperiodic{2});
% 
% index.node_purely_infinite = setdiff(index.node_infinite{1}, ...
%                                        union(index.node_corner{1}, index.node_corner{2}));
% index.element_purely_infinite = setdiff(index.element_infinite{1}, ...
%                                            union(index.element_corner{1}, index.element_corner{2}));
% index.node_purely_interior = setdiff(1 : size(nodes, 2), index.node_all_boundary);
% index.element_purely_interior = setdiff(1 : size(elements, 2), index.element_all_boundary);
% 
% [~, index_sort_1] = sort(nodes(1, index.node_quasiperiodic{1}), 'ascend');
% [~, index_sort_2] = sort(nodes(1, index.node_quasiperiodic{2}), 'ascend');
% index.node_sort_quasiperiodic = cell(2, 1);
% index.node_sort_quasiperiodic{1} = index.node_quasiperiodic{1}(index_sort_1);
% index.node_sort_quasiperiodic{2} = index.node_quasiperiodic{2}(index_sort_2);
% 
% [~, index_sort_infinite] = sort(nodes(2, index.node_infinite{1}), 'ascend');
% [~, index_sort_near_infinite] = sort(nodes(2, index.node_near_infinite{1}), 'ascend');
% [~, index_sort_infinite_opposite] = sort(nodes(2, index.node_infinite_opposite{1}), 'ascend');
% index.node_sort_infinite = cell(1, 1);
% index.node_sort_infinite{1} = index.node_infinite{1}(index_sort_infinite);
% index.node_sort_near_infinite = cell(1, 1);
% index.node_sort_near_infinite{1} = index.node_near_infinite{1}(index_sort_near_infinite);
% index.node_sort_infinite_opposite = cell(1, 1);
% index.node_sort_infinite_opposite{1} = index.node_infinite_opposite{1}(index_sort_infinite_opposite);
% 
% index.node_index_purely_infinite = setdiff(index.node_sort_infinite{1}, node_index_quasiperiodic{2}, 'stable');
% index.node_index_purely_near_infinite = setdiff(index.node_sort_near_infinite{1}, node_index_quasiperiodic{2}, 'stable');
% index.node_index_purely_infinite_opposite = setdiff(index.node_sort_infinite_opposite{1}, node_index_quasiperiodic{2}, 'stable');
% index.node_index_purely_near_infinite_plus_purely_infinite = union(index.node_index_purely_near_infinite, ... 
%                                                                    index.node_index_purely_infinite);
% index.node_index_purely_near_infinite_plus_purely_infinite_opposite = union(index.node_index_purely_near_infinite, ... 
%                                                                    index.node_index_purely_infinite_opposite);
% index.new.free = 1 : size(nodes, 2);
% tmp = union(index.node_index_purely_infinite, index.node_index_purely_near_infinite);
% tmp = union(tmp, index.node_index_purely_infinite_opposite);
% tmp = union(tmp, node_index_quasiperiodic{2});
% index.new.free = setdiff(index.new.free, tmp);
% index.new.free = [index.node_index_purely_infinite_opposite, index.new.free, index.node_index_purely_near_infinite.'];
% % index.new.all = [index.new.free, index.node_index_purely_infinite];
% % index.new.free = setdiff(index.node_index_except_infinite_plus_near_infinite, node_index_quasiperiodic{2});
% % index.new.free = [index.node_index_purely_infinite_opposite, index.new.free];
% index.new.free_end = [index.new.free, index.node_index_purely_infinite];
% 
% mesh.nodes = nodes;
% mesh.elements = elements;
% mesh.elements_simple = elements(1 : 3, :);
% mesh.index = index;
% mesh.MaxElementSize = model.Mesh.MaxElementSize;
% mesh.MinElementSize = model.Mesh.MinElementSize;
% mesh.GeometricOrder = model.Mesh.GeometricOrder;
% 
% mesh.free_nodes = nodes(:, index.new.free);
% % mesh.free_elements = elements(:, index.free_elements);
% mesh.purely_interior_elements = elements(:, index.element_purely_interior);
% mesh.purely_infinite_elements = elements(:, index.element_purely_infinite);
% mesh.purely_infinite_opposite_elements = elements(:, index.element_purely_infinite_opposite);
% mesh.quasiperiodic_elements = cell(2, 1);
% mesh.quasiperiodic_elements{1} = elements(:, index.element_quasiperiodic{1});
% mesh.quasiperiodic_elements{2} = elements(:, index.element_quasiperiodic{2});
% % % find out the locations of boundary nodes in the corresponding elements
% % mesh.index.location_quasiperiodic2_nodes_in_elements = ismember(mesh.quasiperiodic_elements{2}(1 : 3, :), index.node_quasiperiodic{2});
% % mesh.index.location_quasiperiodic2_nodes_in_elements = ismember(mesh.quasiperiodic_elements{2}(1 : 3, :), index.node_quasiperiodic{2});
% 
% mesh.n_nodes = size(nodes, 2);
% mesh.n_elements = size(elements, 2);
% mesh.n_free_nodes = size(mesh.free_nodes, 2);
% % mesh.n_free_elements = size(mesh.free_elements, 2);
% mesh.n_purely_interior_elements = length(index.element_purely_interior);
% mesh.n_purely_infinite_elements = length(index.element_purely_infinite);
% mesh.n_purely_infinite_opposite_elements = length(index.element_purely_infinite_opposite);
% mesh.n_quasiperiodic_elements = zeros(2, 1);
% mesh.n_quasiperiodic_elements(1) = length(index.element_quasiperiodic{1});
% mesh.n_quasiperiodic_elements(2) = length(index.element_quasiperiodic{2});


% eval(['if actual_index_edge_' num2str(i) ' == 1']);
% eval(['    edge_index_for_boundary_1 = actual_index_edge_' num2str(i) ';']);
% eval(['elseif actual_index_edge_' num2str(i) ' == 2']);
% eval(['    edge_index_for_boundary_2 = actual_index_edge_' num2str(i) ';']);
% eval(['elseif actual_index_edge_' num2str(i) ' == 3']);
% eval(['    edge_index_for_boundary_3 = actual_index_edge_' num2str(i) ';']);
% eval(['elseif actual_index_edge_' num2str(i) ' == 4']);
% eval(['    edge_index_for_boundary_4 = actual_index_edge_' num2str(i) ';']);
% eval('end');







% function mesh = classifyVertices(model, parameters)
% 
% geom = model.Geometry;
% nodes = model.Mesh.Nodes;
% elements = model.Mesh.Elements;
% n_boundary = length(parameters.boundary_condition);
% vertex = parameters.geometry.boundary.polygon.vertex;
% 
% % a1 = parameters.a1;
% % a2 = parameters.a2;
% % edge1 = [0, a1(1); 0, a1(2)];
% % edge2 = [a1(1), a1(1) + a2(1); a1(2), a1(2) + a2(2)];
% % edge3 = [a1(1) + a2(1), a2(1); a1(2) + a2(2), a2(2)];
% % edge4 = [a2(1), 0; a2(2), 0];
% edge1 = vertex(1 : 2, :).';
% edge2 = vertex(2 : 3, :).';
% edge3 = vertex(3 : 4, :).';
% edge4 = vertex([4 1], :).';
% 
% %% Edges
% % indexed counterclockwisely, the first edge is 0 -> a1
% index_edge = cell(geom.NumEdges, 1);
% actual_index_edge = zeros(geom.NumEdges, 1);
% for i = 1 : geom.NumEdges
%     index_edge{i} = findNodes(model.Mesh, "region", "Edge", i);
%     actual_index_edge(i) = isOnBoundary(nodes(:, index_edge{i}(2)), edge1, edge2, edge3, edge4);
% end
% edge_index_for_outer_boundary = zeros(n_boundary, 1);
% node_index_for_outer_boundary = cell(n_boundary, 1); 
% for i = 1 : n_boundary
%     edge_index_for_outer_boundary(i) = find(actual_index_edge == i);
%     node_index_for_outer_boundary{i} = index_edge{edge_index_for_outer_boundary(i)};
% end
% edge_index_for_inner_boundary = setdiff(1 : geom.NumEdges, edge_index_for_outer_boundary);
% node_index_for_inner_boundary = cell(geom.NumEdges - n_boundary, 1);
% for i = 1 : geom.NumEdges - n_boundary
%     node_index_for_inner_boundary{i} = index_edge{edge_index_for_inner_boundary(i)};
% end
% 
% % %% relabel indices of nodes and elements
% % flag_inf = 0;
% % for i = 1 : n_boundary
% %     if strcmp(parameters.boundary_condition{i}, 'infinite')
% %         flag_inf = i;
% %     end
% % end
% % if flag_inf ~= 0
% %     node_index_inf = node_index_for_outer_boundary{flag_inf};
% %     element_index_inf = find((sum(ismember(elements, node_index_inf), 1) > 0) ~= 0);
% %     node_index_near_inf = elements(:, element_index_inf);
% %     node_index_near_inf = unique(node_index_near_inf);
% %     node_index_near_inf = setdiff(node_index_near_inf, node_index_inf);
% %     if flag_inf - 2 > 0
% %         flag_inf2 = flag_inf - 2;
% %     else
% %         flag_inf2 = flag_inf + 2;
% %     end
% %     node_index_inf2 = node_index_for_outer_boundary{flag_inf2};
% %     element_index_inf2 = find((sum(ismember(elements, node_index_inf2), 1) > 0) ~= 0);
% %     node_index_near_inf2 = elements(:, element_index_inf2);
% %     node_index_near_inf2 = unique(node_index_near_inf2);
% %     node_index_near_inf2 = setdiff(node_index_near_inf2, node_index_inf2);
% % end
% % 
% % index_all = 1 : size(nodes, 2);
% % index_need_relabel = [node_index_inf2, node_index_near_inf2.', node_index_near_inf.', node_index_inf];
% % index_new = setdiff(index_all, index_need_relabel);
% % index_new = [node_index_inf2, node_index_near_inf2.', index_new, node_index_near_inf.', node_index_inf];
% % [~, index_new] = sort(index_new);
% % nodes = nodes(:, index_new);
% % [~, ~, ic] = unique(elements);
% % elements(:) = index_new(ic);
% 
% % QBC, find nodes near the second boundary with QBC (positive direction)
% flag_QBC = 0;
% for i = 1 : n_boundary
%     if strcmp(parameters.boundary_condition{i}, 'quasi-periodic')
%         flag_QBC = i; 
%     end
%     break; % need only the first one
% end
% if flag_QBC ~= 0
%     node_index_QBC = node_index_for_outer_boundary{flag_QBC};
%     % node_index_near_QBC = [];
%     % for i = 1 : length(node_index_QBC)
%     %     node_index_near_QBC = [node_index_near_QBC; reshape(elements(:, ceil(find(elements == node_index_QBC(i)) / 3)), [], 1)];
%     % end
%     element_index_QBC = find((sum(ismember(elements, node_index_QBC), 1) > 0) ~= 0);
%     node_index_near_QBC = elements(:, element_index_QBC);
%     node_index_near_QBC = unique(node_index_near_QBC);
%     node_index_near_QBC = setdiff(node_index_near_QBC, node_index_QBC);
%     if flag_QBC - 2 > 0
%         flag_QBC2 = flag_QBC - 2;
%     else
%         flag_QBC2 = flag_QBC + 2;
%     end
%     node_index_QBC2 = node_index_for_outer_boundary{flag_QBC2};
%     element_index_QBC2 = find((sum(ismember(elements, node_index_QBC2), 1) > 0) ~= 0);
%     node_index_near_QBC2 = elements(:, element_index_QBC2);
%     node_index_near_QBC2 = unique(node_index_near_QBC2);
%     node_index_near_QBC2 = setdiff(node_index_near_QBC2, node_index_QBC2);
% end
% 
% % Neumann, the same operation as that for QBC
% flag_Neumann = [];
% for i = 1 : n_boundary
%     if strcmp(parameters.boundary_condition{i}, 'Neumann')
%         flag_Neumann = [flag_Neumann; i];
%     end
% end
% node_index_Neumann = cell(length(flag_Neumann), 1);
% node_index_near_Neumann = cell(length(flag_Neumann), 1);
% element_index_Neumann = cell(length(flag_Neumann), 1);
% for i = 1 : length(flag_Neumann)
%     node_index_Neumann{i} = node_index_for_outer_boundary{flag_Neumann(i)};
%     element_index_Neumann{i} = find((sum(ismember(elements, node_index_Neumann{i}), 1) > 0) ~= 0);
%     node_index_near_Neumann{i} = elements(:, element_index_Neumann{i});
%     node_index_near_Neumann{i} = unique(node_index_near_Neumann{i});
%     node_index_near_Neumann{i} = setdiff(node_index_near_Neumann{i}, node_index_Neumann{i});
% end
% 
% % Dirichlet, the same operation as that for QBC
% flag_Dirichlet = [];
% for i = 1 : n_boundary
%     if strcmp(parameters.boundary_condition{i}, 'Dirichlet')
%         flag_Dirichlet = [flag_Neumann; i];
%     end
% end
% node_index_Dirichlet = cell(length(flag_Dirichlet), 1);
% node_index_near_Dirichlet = cell(length(flag_Dirichlet), 1);
% element_index_Dirichlet = cell(length(flag_Dirichlet), 1);
% for i = 1 : length(flag_Dirichlet)
%     node_index_Dirichlet{i} = node_index_for_outer_boundary{flag_Dirichlet(i)};
%     element_index_Dirichlet{i} = find((sum(ismember(elements, node_index_Dirichlet{i}), 1) > 0) ~= 0);
%     node_index_near_Dirichlet{i} = elements(:, element_index_Dirichlet{i});
%     node_index_near_Dirichlet{i} = unique(node_index_near_Dirichlet{i});
%     node_index_near_Dirichlet{i} = setdiff(node_index_near_Dirichlet{i}, node_index_Dirichlet{i});
% end
% 
% % infinite, the same operation as that for QBC
% flag_inf = 0;
% for i = 1 : n_boundary
%     if strcmp(parameters.boundary_condition{i}, 'infinite')
%         flag_inf = i;
%     end
% end
% if flag_inf ~= 0
%     node_index_inf = node_index_for_outer_boundary{flag_inf};
%     % node_index_near_inf = [];
%     % for i = 1 : length(node_index_inf)
%     %     node_index_near_inf = [node_index_near_inf; reshape(elements(:, ceil(find(elements == node_index_inf(i)) / 3)), [], 1)];
%     % end
%     element_index_inf = find((sum(ismember(elements, node_index_inf), 1) > 0) ~= 0);
%     node_index_near_inf = elements(:, element_index_inf);
%     node_index_near_inf = unique(node_index_near_inf);
%     node_index_near_inf = setdiff(node_index_near_inf, node_index_inf);
%     if flag_inf - 2 > 0
%         flag_inf2 = flag_inf - 2;
%     else
%         flag_inf2 = flag_inf + 2;
%     end
%     node_index_inf2 = node_index_for_outer_boundary{flag_inf2};
%     % node_index_near_inf2 = [];
%     % for i = 1 : length(node_index_inf2)
%     %     node_index_near_inf2 = [node_index_near_inf2; reshape(elements(:, ceil(find(elements == node_index_inf2(i)) / 3)), [], 1)];
%     % end
%     element_index_inf2 = find((sum(ismember(elements, node_index_inf2), 1) > 0) ~= 0);
%     node_index_near_inf2 = elements(:, element_index_inf2);
%     node_index_near_inf2 = unique(node_index_near_inf2);
%     node_index_near_inf2 = setdiff(node_index_near_inf2, node_index_inf2);
% end
% 
% 
% %% Elements
% % the first index means the indices of nodes and elements located in the air
% index_1 = find(nodes(1, :) == 0);
% index_2 = find(nodes(2, :) == 0);
% index_0 = intersect(index_1, index_2);
% node_index_for_face = cell(geom.NumFaces, 1);
% flag_face = 0;
% for i = 1 : geom.NumFaces
%     node_index_for_face{i} = findNodes(model.Mesh, "region", "Face", i);
%     if ismember(index_0, node_index_for_face{i})
%         flag_face = i;
%     end
% end
% if flag_face ~= 1
%     tmp = node_index_for_face{1};
%     node_index_for_face{1} = node_index_for_face{flag_face};
%     node_index_for_face{flag_face} = tmp;
% end
% 
% elements = [elements; zeros(1, size(elements, 2))]; % the fourth row stores the indices of faces where the element lies in
% element_index_for_face = cell(geom.NumFaces, 1);
% flag_face = 0;
% for i = 1 : geom.NumFaces
%     element_index_for_face{i} = findElements(model.Mesh, "region", "Face", i);
%     elements(4, element_index_for_face{i}) = i;
%     if ismember(index_0, elements(1 : 3, element_index_for_face{i}))
%         flag_face = i;
%     end
% end
% if flag_face ~= 1
%     tmp = element_index_for_face{1};
%     element_index_for_face{1} = node_index_for_face{flag_face};
%     element_index_for_face{flag_face} = tmp;
% 
%     elements(4, element_index_for_face{1}) = flag_face;
%     elements(4, element_index_for_face{flag_face}) = 1;
% end
% 
% index.node_for_inner_boundary = node_index_for_inner_boundary;
% index.node_for_outer_boundary = node_index_for_outer_boundary;
% index.free_nodes = 1 : size(nodes, 2);
% index.free_elements = 1 : size(elements, 2);
% if ~isempty(flag_QBC)
%     index.node_QBC = node_index_QBC;
%     index.node_near_QBC = node_index_near_QBC;
%     index.element_QBC = element_index_QBC;
%     index.node_QBC2 = node_index_QBC2;
%     index.node_near_QBC2 = node_index_near_QBC2;
%     index.element_QBC2 = element_index_QBC2;
%     % free nodes do not include the second QBC
%     index.free_nodes = setdiff(index.free_nodes, node_index_QBC2);
%     index.free_elements = setdiff(index.free_elements, element_index_QBC2);
% end
% if ~isempty(flag_Neumann)
%     index.node_Neumann = node_index_Neumann;
%     index.node_near_Neumann = node_index_near_Neumann;
%     index.element_Neumann = element_index_Neumann;
% end
% if ~isempty(flag_Dirichlet)
%     index.node_Dirichlet = node_index_Dirichlet;
%     index.node_near_Dirichlet = node_index_near_Dirichlet;
%     index.element_Dirichlet = element_index_Dirichlet;
% end
% if flag_inf ~= 0
%     index.node_inf = node_index_inf;
%     index.node_inf2 = node_index_inf2;
%     index.node_near_inf = node_index_near_inf;
%     index.node_near_inf2 = node_index_near_inf2;
%     index.element_inf = element_index_inf;
%     index.element_inf2 = element_index_inf2;
%     % free nodes do not include the first inf (actual inf)
%     index.free_nodes = setdiff(index.free_nodes, node_index_inf);
%     index.free_elements = setdiff(index.free_elements, element_index_inf);
% end
% index.node_for_face = node_index_for_face;
% index.element_for_face = element_index_for_face;
% index.element_pure_QBC2 = setdiff(index.element_QBC2, index.element_inf); % OK
% index.element_pure_inf = setdiff(index.element_inf, index.element_QBC2); % OK
% index.element_pure_inf2 = setdiff(index.element_inf2, index.element_QBC2); % OK
% index.element_pure_interior = setdiff(1 : size(elements, 2), ... 
%     union(union(index.element_pure_inf, index.element_pure_inf2), index.element_QBC2));
% 
% mesh.nodes = nodes;
% mesh.elements = elements;
% mesh.index = index;
% mesh.MaxElementSize = model.Mesh.MaxElementSize;
% mesh.MinElementSize = model.Mesh.MinElementSize;
% mesh.GeometricOrder = model.Mesh.GeometricOrder;
% 
% mesh.free_nodes = nodes(:, index.free_nodes);
% mesh.free_elements = elements(:, index.free_elements);
% mesh.QBC_nodes = nodes(:, index.node_QBC);
% mesh.QBC_elements = elements(:, index.element_QBC);
% mesh.QBC2_nodes = nodes(:, index.node_QBC2);
% mesh.QBC2_elements = elements(:, index.element_QBC2);
% mesh.pure_QBC2_elements = elements(:, index.element_pure_QBC2);
% mesh.inf_nodes = nodes(:, index.node_inf);
% mesh.near_inf_nodes = nodes(:, index.node_near_inf);
% mesh.inf_elements = elements(:, index.element_inf);
% mesh.pure_inf_elements = elements(:, index.element_pure_inf);
% mesh.pure_inf2_elements = elements(:, index.element_pure_inf2);
% mesh.pure_interior_elements = elements(:, index.element_pure_interior);
% 
% mesh.n_nodes = size(nodes, 2);
% mesh.n_free_nodes = length(index.free_nodes);
% mesh.n_elements = size(elements, 2);
% mesh.n_free_elements = length(index.free_elements);
% mesh.n_QBC_nodes = length(index.node_QBC);
% mesh.n_QBC_elements = length(index.element_QBC);
% mesh.n_QBC2_nodes = length(index.node_QBC2);
% mesh.n_QBC2_elements = length(index.element_QBC2);
% mesh.n_pure_QBC2_elements = length(index.element_pure_QBC2);
% mesh.n_inf_nodes = length(index.node_inf);
% mesh.n_inf2_nodes = length(index.node_inf2);
% mesh.n_near_inf_nodes = length(index.node_near_inf);
% mesh.n_near_inf2_nodes = length(index.node_near_inf2);
% mesh.n_inf_elements = length(index.element_inf);
% mesh.n_inf2_elements = length(index.element_inf2);
% mesh.n_pure_inf_elements = length(index.element_pure_inf);
% mesh.n_pure_inf2_elements = length(index.element_pure_inf2);
% mesh.n_pure_interior_elements = length(index.element_pure_interior);
% 
% %% sort indices with 
% % 1. nodes with QBC to be abandoned
% % 2. free nodes
% % 3. nodes with inf BC
% index_all = 1 : size(nodes, 2);
% index_need_relabel = unique([node_index_near_inf.', node_index_inf], 'stable');
% index_new = setdiff(index_all, index_need_relabel);
% index_QBC = setdiff(node_index_QBC, node_index_inf);
% index_QBC2 = setdiff(node_index_QBC2, node_index_inf);
% index_note = intersect(node_index_QBC2, node_index_near_inf.');
% index_tmp = setdiff(node_index_QBC, index_need_relabel);
% index_tmp2 = setdiff(node_index_QBC2, index_need_relabel);
% index_free = setdiff(index_new, union(index_tmp, index_tmp2));
% index_QBC2_and_inf = intersect(node_index_QBC2, node_index_inf);
% index_new = [index_tmp2, index_tmp, index_free, setdiff(index_need_relabel, index_QBC2_and_inf)];
% 
% new_index_inf2 = setdiff(node_index_inf2, intersect(node_index_inf2, node_index_QBC2), 'stable');
% new_index_inf = setdiff(node_index_inf, intersect(node_index_inf, node_index_QBC2), 'stable');
% new_index_near_inf = setdiff(node_index_near_inf.', intersect(node_index_near_inf.', node_index_QBC2), 'stable');
% index_else = union(new_index_inf, new_index_near_inf);
% index_else = union(index_else, new_index_inf2);
% index_else = union(index_else, node_index_QBC2);
% index_free = setdiff(index_all, index_else);
% % make sure the order is right
% index_free = [new_index_inf2, index_free, new_index_near_inf];
% index_free_end = [index_free, new_index_inf];
% mesh.new.new_index_inf = new_index_inf;
% mesh.new.new_index_inf2 = new_index_inf2;
% mesh.new.new_index_near_inf = new_index_near_inf;
% 
% % [~, index_new2] = sort(index_new);
% % nodes_new = nodes(:, index_new2);
% % [~, ~, ic] = unique(elements(1 : 3, :));
% % elements_new = reshape(index_new2(ic), 3, []);
% % mesh.new.nodes = nodes_new;
% % mesh.new.elements = elements_new;
% mesh.new.index = index_new;
% % mesh.new.index2 = index_new2;
% mesh.new.index_free = index_free;
% mesh.new.index_free_end = index_free_end;
% mesh.n_end_free_nodes = length(index_free_end);
% % reorder index_QBC2 corresponding to index_QBC
% new_node_QBC = nodes(:, index_QBC);
% new_node_QBC2 = nodes(:, index_QBC2);
% [~, index_sort] = sort(new_node_QBC(1, :), 'ascend');
% [~, index_sort2] = sort(new_node_QBC2(1, :), 'ascend');
% mesh.new.index_QBC = index_QBC(index_sort);
% mesh.new.index_QBC2 = index_QBC2(index_sort2);
% mesh.new.index_QBC2_and_inf = index_QBC2_and_inf;
% mesh.new.index_inf = node_index_inf;
% mesh.new.index_near_QBC2 = setdiff(node_index_near_QBC2.', node_index_inf);
% mesh.new.index_QBC_plus_near_QBC2 = [mesh.new.index_QBC, mesh.new.index_near_QBC2];
% mesh.new.index_QBC2_plus_near_QBC2 = [mesh.new.index_QBC2, mesh.new.index_near_QBC2];
% 
% mesh.new.index_inf = node_index_inf;
% mesh.new.index_inf = setdiff(mesh.new.index_inf, node_index_QBC2);
% mesh.new.index_inf2 = node_index_inf2;
% mesh.new.index_inf2 = setdiff(mesh.new.index_inf2, node_index_QBC2);
% mesh.new.index_near_inf = setdiff(node_index_near_inf.', node_index_QBC2, 'stable');
% mesh.new.index_inf_plus_near_inf = [new_index_inf, new_index_near_inf];
% mesh.new.index_near_inf_plus_inf = [new_index_near_inf, new_index_inf];
% new_node_QBC_all = nodes(:, node_index_QBC);
% new_node_QBC2_all = nodes(:, node_index_QBC2);
% [~, index_sort_all] = sort(new_node_QBC_all(1, :), 'ascend');
% [~, index_sort2_all] = sort(new_node_QBC2_all(1, :), 'ascend');
% mesh.new.index_QBC_all = node_index_QBC(index_sort_all);
% mesh.new.index_QBC2_all = node_index_QBC2(index_sort2_all);
% mesh.new.index_QBC_plus_near_QBC2_all = [mesh.new.index_QBC_all, node_index_near_QBC2.'];
% mesh.new.index_QBC2_plus_near_QBC2_all = [mesh.new.index_QBC2_all, node_index_near_QBC2.'];
% 
% 
% mesh.flag.QBC = flag_QBC;
% mesh.flag.QBC2 = flag_QBC2;
% mesh.flag.Neumann = flag_Neumann;
% mesh.flag.inf = flag_inf;
% mesh.flag.inf2 = flag_inf2;
% 
% end 
% 
% 
% 
% 
% 
% 
% % eval(['if actual_index_edge_' num2str(i) ' == 1']);
% % eval(['    edge_index_for_boundary_1 = actual_index_edge_' num2str(i) ';']);
% % eval(['elseif actual_index_edge_' num2str(i) ' == 2']);
% % eval(['    edge_index_for_boundary_2 = actual_index_edge_' num2str(i) ';']);
% % eval(['elseif actual_index_edge_' num2str(i) ' == 3']);
% % eval(['    edge_index_for_boundary_3 = actual_index_edge_' num2str(i) ';']);
% % eval(['elseif actual_index_edge_' num2str(i) ' == 4']);
% % eval(['    edge_index_for_boundary_4 = actual_index_edge_' num2str(i) ';']);
% % eval('end');



end