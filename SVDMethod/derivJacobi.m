function dcdx=derivJacobi(state,mu)

%C=jacobiConst(states,mu)
%
% Computes the jacobi constant for a given states
%
% inputs:
%   states: (x,y,z=0, xv,yv,zv)
%   mu: mass parameter
%
% outputs:
%   dcdx: jacobi derivative

x = state(1);
y = state(2);
z = state(3);
xv = state(4);
yv = state(5);
zv = state(6);

%the distances
r1=sqrt( (x+mu)^2 + y^2 + z^2 );
r2=sqrt( (x-(1-mu))^2 + y^2 + z^2);

%Compute the Jacobi Energy
% C=-(xv^2 + yv^2 + zv^2)/2 + 2*( (x^2 + y^2)/2 + (1-mu)/r1 + mu/r2);

Cx = x - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3;
Cy = y - (1-mu)*y/r1^3 - mu*y/r2^3;
Cxv = -xv;
Cyv = -yv;
Czv = -zv;

dcdx = [Cx Cy Cxv Cyv Czv];

end