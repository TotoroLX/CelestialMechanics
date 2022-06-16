function d=missDist(x0, x1, C0, C1)

%d=missDis(x0, x1, C0, C1)
%
% Computes normalized miss distance
%
% inputs:
%   x0: initial state
%   x1: final state
%   C0: initial jacobi constant
%   C1: final jacobi constatn
%
% outputs:
%   d: normalized diss distance

v0 = x0(4:6);
v1 = x1(4:6);

x0 = x0(1:3);
x1 = x1(1:3);

d = norm(x1-x0)/norm(x0) + norm(v1-v0)/norm(v0) + norm(C1/C0-1);

end