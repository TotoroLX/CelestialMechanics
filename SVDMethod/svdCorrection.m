function dKdxi0Inv = svdCorrection(phi, xT, mu, dcdx, epsilon)

%dx0 = svdCorrection(phi, xT, mu)
%
% SVD correction by Ryan P. Russell and M Lara
%
% inputs:
%   phi: state transition matrix
%   xT: state at t = T
%   mu: mass paramter
%   dcdx: Jacobi derivative
%   epsilon: parameter
%
% outputs:
%   dKdxi0Inv: variation of initial state

global param

param = mu;

dxT = cr3bp_svd(0, xT);
dzT = dxT(3);

idx = [1 2 4 5 6];

dxiT = dxT(idx);

phi_3 = phi(3,idx);
phi_r = phi(idx, idx);

dxTdx0 = phi_r - 1/dzT*dxiT*phi_3;

dKdxi0 = [dxTdx0 - eye(5)];%; dcdx];

[U,D,V] = svd(dKdxi0);

S = zeros(size(D'));
for i=1:5
    if D(i,i) <= epsilon
        S(i,i) = 0;
    else
        S(i,i) = 1/D(i,i);
    end
end

dKdxi0Inv = V*S*U';

end
