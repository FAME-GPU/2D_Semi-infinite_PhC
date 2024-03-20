function [S, detJ] = isoparametricMapVectorize(x, y, r, s, shapefcn)

[S, dSdr, dSds] = shapefcn(r, s);
j11 = dSdr.' * x;  % \partial x / \partial r
j12 = dSdr.' * y;  % \partial y / \partial r
j21 = dSds.' * x;  % \partial x / \partial s
j22 = dSds.' * y;  % \partial y / \partial s
detJ = j11 .* j22 - j12 .* j21;

% dSdx = ( j22 .* dSdr - j12 .* dSds) ./ detJ;
% dSdy = (-j21 .* dSdr + j11 .* dSds) ./ detJ;

end