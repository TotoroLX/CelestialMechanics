function PHIdot = vareqn3bp(t, PHI)

%PHIdot = vareqn3bp(t, PHI)
%
% Differential equation for state transition (or sensitivity) matrix 
%
% d PHI(t,t0)
% ----------- = DU(t) * PHI(t,t0)
%      dt
%
% See Ch 5.12 and Eq 5.190 By Prof. Casotto

global param

mu = param;

m1 = 1-mu; % LARGER MASS
m2 = mu; % smaller mass
N = 6;

x(1:N) = PHI(1+N^2:N+N^2);
%phi = reshape(PHI(1:N^2),N,N);
%Make phi
for i=1:6
    for j=1:6
        phi(i,j) = PHI(6*(i-1)+j);
    end
end

r1 = sqrt( (PHI(37)+mu)^2 + PHI(38)^2 + PHI(39)^2 ); % r1 distance to m1
r2 = sqrt( (PHI(37)-(1-mu))^2 + PHI(38)^2 + PHI(39)^2 ); % r2 distance to m2

Ux = PHI(37) - (1-mu)*(PHI(37)+mu)/r1^3 - mu*(PHI(37)-(1-mu))/r2^3;
Uy = PHI(38) - (1-mu)*PHI(38)/r1^3 - mu*PHI(38)/r2^3;
Uz = -(1-mu)*PHI(39)/r1^3 - mu*PHI(39)/r2^3;

m1p = -mu; % position of LARGER MASS
m2p = 1-mu; % position of smaller mass

Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 +...
    3*(1-mu)*(PHI(37)-m1p)^2/r1^5 + 3*mu*(PHI(37)-m2p)^2/r2^5;
Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*((1-mu)/r1^5 + mu/r2^5)*PHI(38)^2;
Uzz = - (1-mu)/r1^3 - mu/r2^3 + 3*((1-mu)/r1^5 + mu/r2^5)*PHI(39)^2;
Uxy = 3*((1-mu)*(PHI(37)-m1p)/r1^5 + mu*(PHI(37)-m2p)/r2^5)*PHI(38);
Uxz = 3*((1-mu)*(PHI(37)-m1p)/r1^5 + mu*(PHI(37)-m2p)/r2^5)*PHI(39);
Uyz = 3*((1-mu)/r1^5 + mu/r2^5)*PHI(38)*PHI(39);

DU = [0   0   0   1  0 0;...
      0   0   0   0  1 0;...
      0   0   0   0  0 1;...
      Uxx Uxy Uxz 0  2 0;...
      Uxy Uyy Uyz -2 0 0;...
      Uxz Uyz Uzz 0  0 0];

phidot = DU * phi;

PHIdot = zeros(N^2+N,1);
for i=1:6
    for j=1:6
        PHIdot(6*(i-1)+j) = phidot(i,j);
    end
end
PHIdot(N^2+1) = PHI(40);
PHIdot(N^2+2) = PHI(41);
PHIdot(N^2+3) = PHI(42);
PHIdot(N^2+4) = 2*PHI(41) + Ux;
PHIdot(N^2+5) = -2*PHI(40) + Uy;
PHIdot(N^2+6) = Uz;

end