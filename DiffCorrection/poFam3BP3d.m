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

x0poGuess1=[1.176999999999998;0;0.059682484925404;0;-0.173565964236403;0];
% x0=[0.823389831024856;0;0.003724454349521;0;0.126557411688950;0];
x0poGuess2=[1.177999999999998;0;0.055682484925404;0;-0.170565964236403;0];
% x0=[1.000000000000000;0;0.056667402379709;0;0.591715153162376;0]; 2*1.310229401292405;
% x0=[1.010000000000000;0;0.065589811196310;0;0.523743735947168;0]; 2*1.403831515057635;
x0poGuess1=[-0.160000000000000; 0; 1.913211277268421; 0; 0.117265264191552; 0];
x0poGuess2=[-0.150000000000000; 0; 1.909867588631233; 0; 0.109735057876700; 0];

% First guess
iFam = 1 ;
fprintf('::poFamGet : number %d\n',iFam) ;
%[x01, t1] = poDifCor3BP3d(x0poGuess1, mu) ;
[x01, t1] = poDifCor3BP3d(x0poGuess1, mu) ;

% Second guess
iFam = 2;
fprintf('::poFamGet : number %d\n',iFam) ;
%[x02, t2] = poDifCor3BP3d(x0poGuess2, mu) ;
[x02, t2] = poDifCor3BP3d(x0poGuess2, mu) ;

x0po(1:2, 1:N) = [x01(:)'; x02(:)'];
T(1:2) = [2*t1; 2*t2];

for iFam = 3:nFam
    fprintf('::poFamGet : number %d\n',iFam);
    
    dx = x0po(iFam-1,1) - x0po(iFam-2, 1);
    dz = x0po(iFam-1,3) - x0po(iFam-2, 3);
    dzp = x0po(iFam-1,5) - x0po(iFam-2, 5);
    
    dt = T(iFam-1) - T(iFam-2);
    
    % new initial guess
    x0po_g = [x0po(iFam-1,1) + dx; 0;...
              x0po(iFam-1,3) + dz; 0; x0po(iFam-1,5) + dzp; 0];
    t1po_g = (T(iFam-1) + dt)/2 + delt;
    
    % differential correction
    %[x0po_iFam, t1_iFam] = poDifCor3BP3d(x0po_g, mu);%, t1po_g) ;
    [x0po_iFam, t1_iFam] = poDifCor3BP3d(x0po_g, mu);%, t1po_g) ;
    
    x0po(iFam, 1:N) = x0po_iFam' ;
    T(iFam) = 2*t1_iFam ;
end

% Store the family of initial states
dum = [x0po T];
% save x0po_T.dat -ascii -double dum

end