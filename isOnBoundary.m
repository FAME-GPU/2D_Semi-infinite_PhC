function is_on_boundary = isOnBoundary(node, edge1, edge2, edge3, edge4)

is_on_boundary = 0;

tmp1 = [edge1 - node; 0, 0];
tmp2 = [edge2 - node; 0, 0];
tmp3 = [edge3 - node; 0, 0];
tmp4 = [edge4 - node; 0, 0];
if norm(cross(tmp1(:, 1), tmp1(:, 2))) < 1e-14
    is_on_boundary = 1;
end
if norm(cross(tmp2(:, 1), tmp2(:, 2))) < 1e-14
    is_on_boundary = 2;
end
if norm(cross(tmp3(:, 1), tmp3(:, 2))) < 1e-14
    is_on_boundary = 3;
end
if norm(cross(tmp4(:, 1), tmp4(:, 2))) < 1e-14
    is_on_boundary = 4;
end


end




% function is_on_boundary = isOnBoundary(nodes, a1, a2)
% 
% n = size(nodes, 2);
% if size(nodes, 2) <= 2
%     error('Two few vertices.');
% end
% 
% edge1 = [0, a1(1); 0, a1(2)];
% edge2 = [0, a2(1); 0, a2(2)];
% edge3 = [a1(1), a1(1) + a2(1); a1(2), a1(2) + a2(2)];
% edge4 = [a2(1), a1(1) + a2(1); a2(2), a1(2) + a2(2)];
% is_on_boundary = zeros(n - 2, 1);
% 
% for i = 2 : n - 1
%     tmp1 = [edge1 - nodes(:, i); 0, 0];
%     tmp2 = [edge2 - nodes(:, i); 0, 0];
%     tmp3 = [edge3 - nodes(:, i); 0, 0];
%     tmp4 = [edge4 - nodes(:, i); 0, 0];
%     if norm(cross(tmp1(:, 1), tmp1(:, 2))) < 1e-14
%         is_on_boundary(i - 1) = 1;
%     end
%     if norm(cross(tmp2(:, 1), tmp2(:, 2))) < 1e-14
%         is_on_boundary(i - 1) = 2;
%     end
%     if norm(cross(tmp3(:, 1), tmp3(:, 2))) < 1e-14
%         is_on_boundary(i - 1) = 3;
%     end
%     if norm(cross(tmp4(:, 1), tmp4(:, 2))) < 1e-14
%         is_on_boundary(i - 1) = 4;
%     end
% end
% 
% if norm(is_on_boundary - is_on_boundary(1), inf) == 0
%     is_on_boundary = is_on_boundary(1);
% else
%     is_on_boundary = 0;
% end
% 
% end
% 
% 
