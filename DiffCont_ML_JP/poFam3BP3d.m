function [x0po, T] = poFam3BP3d(mu, eqNum, Ax1, Ax2, nFam)

%[x0po, T] = poFam3BP3d(mu, eqNum, Ax1, Ax2, nFam)
% 
% Generates the initial states for a family of periodic orbits
%
% inputs:
%   mu: mass parameter
%   eqNum: the number of equilibrium point
%   Ax1, Ax2: used for initial guess
%   nFam: number of periodic orbits
%
% outputs:
%   x0po: family of initial states of periodic orbits
%   T: correspond periods

% delt = guessed change in period between successive orbits in family
delt = -1.e-2 ; % <==== may need to be changed

N = 6; % dimension of phase space

x0po = zeros(nFam, N);
T = zeros(nFam, 1);

[x0poGuess1, TGuess1] = poGuessLinear3BP3d(mu, eqNum, Ax1);
[x0poGuess2, TGuess2] = poGuessLinear3BP3d(mu, eqNum, Ax2);
x0poGuess2=[0.989290839334460;0;0.014;0;-1.293979714378937;0];
x0poGuess1=[0.990290839334460;0;0.015;0;-1.303979714378937;0];

iFam = 1 ;
fprintf('::poFamGet : number %d\n',iFam) ;
[x01, t1] = poDifCor3BP3d(x0poGuess1, mu) ;
% x01 = [0.990179228047932;0;-0.023000000000000;0;0.999572391984704;0];
% t1 = 1.893992948457965/2;
myplot(x01,t1*2,mu);

iFam = 2;
fprintf('::poFamGet : number %d\n',iFam) ;
[x02, t2] = poDifCor3BP3d(x0poGuess2, mu) ;
% x02 = [0.990103769097997;0;-0.022000000000000;0;1.023036548525750;0];
% t2 = 1.883095724403974/2;
myplot(x02,t2*2,mu);

x0po(1:2, 1:N) = [x01(:)'; x02(:)'];
T(1:2) = [2*t1; 2*t2];

for iFam = 3:nFam
    fprintf('::poFamGet : number %d\n',iFam);
    
    dx = x0po(iFam-1,1) - x0po(iFam-2, 1);
    dz = x0po(iFam-1,3) - x0po(iFam-2, 3);
    dyp = x0po(iFam-1,5) - x0po(iFam-2, 5);
    
    dt = T(iFam-1) - T(iFam-2);
    
    % new initial guess
    x0po_g = [x0po(iFam-1,1) + dx*0; 0;         x0po(iFam-1,3) + dz;...
              0;                      x0po(iFam-1,5) + dyp;    0];
    t1po_g = (T(iFam-1) + dt)/2 + delt;
    
    % differential correction
    [x0po_iFam, t1_iFam] = poDifCor3BP3d(x0po_g, mu, t1po_g) ;
    myplot(x0po_iFam,t1_iFam*2,mu);
    
    x0po(iFam, 1:N) = x0po_iFam' ;
    T(iFam) = 2*t1_iFam ;
end

% Store the family of initial states
dum = [x0po T];
save x0po_T_m2_13.dat -ascii -double dum

end