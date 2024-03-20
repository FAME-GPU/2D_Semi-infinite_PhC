function [S, dSdr, dSds] = P1(r, s)

S = reshape([1 - r - s; r; s], length(r), []).';
% dSdr = [-1; 1; 0];
% dSds = [-1; 0; 1];
dSdr = reshape(repmat([-1; 1; 0], length(r), 1), [], length(r));
dSds = reshape(repmat([-1; 0; 1], length(r), 1), [], length(r));

% S2 = [1 - r - s; r; s];
% dSdr2 = [-1; 1; 0];
% dSds2 = [-1; 0; 1];

end