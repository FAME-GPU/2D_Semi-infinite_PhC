function [area, a, b, c] = hatGradients(xx, yy)

area = polyarea(xx, yy);

a = [xx(2, :) .* yy(3, :) - xx(3, :) .* yy(2, :); xx(3, :) .* yy(1, :) - xx(1, :) .* yy(3, :); ... 
     xx(1, :) .* yy(2, :) - xx(2, :) .* yy(1, :)] ./ area;
a = a / 2;
b = [yy(2, :) - yy(3, :); yy(3, :) - yy(1, :); yy(1, :) - yy(2, :)] ./ area;
b = b / 2;
c = [xx(3, :) - xx(2, :); xx(1, :) - xx(3, :); xx(2, :) - xx(1, :)] ./ area;
c = c / 2;

end