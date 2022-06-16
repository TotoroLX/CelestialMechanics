function DU = DUmatrix3BP(x, mu)

%DU = DUmatrix3BP(x, mu)
%
% Computes the Hessian matrix of potential U, where xdot = U_x
% for a given point x in 3-d synodic space
%
% inputs:
%   x: x(1),x(2),x(3) coordinates
%   mu: mass parameter
%
% outputs:
%   DU: Jacobian matrix
%
%-----------------------------------------------------------
% CR3BP CONVENTION:
%         L4
%
%
% L3-----M1-------L1---M2---L2     M1=1-mu, M2=mu
%       -mu           1-mu
%
%         L5
% 
% Follow Ch 5.9 Prof. Casotto

m1p = -mu;
m2p = 1-mu;

r1 = sqrt( (x(1)-m1p)^2 + x(2)^2 + x(3)^2); % distance to M1, LARGER MASS
r2 = sqrt( (x(1)-m2p)^2 + x(2)^2 + x(3)^2); % distance to M2, smaller mass

Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 +...
    3*(1-mu)*(x(1)-m1p)^2/r1^5 + 3*mu*(x(1)-m2p)^2/r2^5;
Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*((1-mu)/r1^5 + mu/r2^5)*x(2)^2;
Uzz = - (1-mu)/r1^3 - mu/r2^3 + 3*((1-mu)/r1^5 + mu/r2^5)*x(3)^2;
Uxy = 3*((1-mu)*(x(1)-m1p)/r1^5 + mu*(x(1)-m2p)/r2^5)*x(2);
Uxz = 3*((1-mu)*(x(1)-m1p)/r1^5 + mu*(x(1)-m2p)/r2^5)*x(3);
Uyz = 3*((1-mu)/r1^5 + mu/r2^5)*x(2)*x(3);

DU = [0   0   0   1  0 0;...
      0   0   0   0  1 0;...
      0   0   0   0  0 1;...
      Uxx Uxy Uxz 0  2 0;...
      Uxy Uyy Uyz -2 0 0;...
      Uxz Uyz Uzz 0  0 0];

end