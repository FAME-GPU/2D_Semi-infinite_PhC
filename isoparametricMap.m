function [S, dSdx, dSdy, detJ] = isoparametricMap(x, y, r, s, shapefcn)

[S, dSdr, dSds] = shapefcn(r, s);
j11 = dot(dSdr, x); 
j12 = dot(dSdr, y);
j21 = dot(dSds, x); 
j22 = dot(dSds, y);
detJ = j11 * j22 - j12 * j21;
dSdx = ( j22 * dSdr - j12 * dSds) / detJ;
dSdy = (-j21 * dSdr + j11 * dSds) / detJ;

end