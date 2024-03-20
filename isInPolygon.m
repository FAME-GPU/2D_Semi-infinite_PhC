function [in, out] = isInPolygon(nodes, vertex)

tol = 1e-5;

v1 = vertex(1, :)';
v2 = vertex(2, :)';
v4 = vertex(4, :)';

h1 = abs(det([v2 - v1, v4 - v1])) / norm(v2 - v1);
h2 = abs(det([v4 - v1, v2 - v1])) / norm(v4 - v1);

dist = zeros(4, size(nodes, 2));
for j = 1 : 4
    v1 = vertex(j, :)';
    v2 = vertex(mod(j, 4) + 1, :)';
    v2mv1 = v2 - v1;

    norm_v2mv1 = norm(v2mv1);
    nodes_minus_v1 = nodes - repmat(v1, 1, size(nodes, 2));
    dist(j, :) = abs(v2mv1(1) * nodes_minus_v1(2, :) - v2mv1(2) * nodes_minus_v1(1, :)) / norm_v2mv1;
end

index_in_1 = dist(1, :) < h1 + tol;
index_in_2 = dist(2, :) < h2 + tol;
index_in_3 = dist(3, :) < h1 + tol;
index_in_4 = dist(4, :) < h2 + tol;
index_in = [index_in_1; index_in_2; index_in_3; index_in_4];
in = find(sum(index_in) == 4);
out = find(sum(index_in) < 4);


end