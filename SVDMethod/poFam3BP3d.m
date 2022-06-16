function [x0po, T] = poFam3BP3d(mu, eqNum, Ax1, Ax2, nFam, epsilon)

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

[x0poGuess1, TGuess1] = poGuessLinear3BP3d(mu, eqNum, Ax2);
[x0poGuess2, TGuess2] = poGuessLinear3BP3d(mu, eqNum, Ax1);
% x0poGuess1 = [0.987080533067444;-0.027216335836602;0;...
%    -0.113218729948413;-0.621170469770948;0.664658830988097]; %DAMP=0.1
% x0poGuess = [0.844641367268382;-0.056452368910337;0;...
%     -0.035700546293609;-0.014036796624513;0.021171332709162]; good guess
x0poGuess1 = [1.399;-0.022;0;...
    0.025022945507023;-0.089351053939445;-0.000022783822845];
x0poGuess1 = [0.987856520000000;-0.003785823100000;0;...
    -0.068141039000000;1.756531000000000;-1.790220000000000];
x0poGuess1 = [0.800000000000000;0;0.574135940790300;...
    0;0.195502355325772;0];


iFam = 1 ;
fprintf('::poFamGet : number %d\n',iFam) ;
[x01, t1] = poSVD3BP3d(x0poGuess1, mu, epsilon);
% x01 = [0.823386384377028;0;-0.008;0;0.127386418509082;0];
% t1 = 1.371716437600570;

x0po(1, 1:N) = x01(:)';
T(1) = t1;

for iFam = 2:nFam
    fprintf('::poFamGet : number %d\n',iFam);
    
%     dx = x0po(iFam-1,1) - x0po(iFam-2, 1);
%     dy = x0po(iFam-1,2) - x0po(iFam-2, 2);
%     dxp = x0po(iFam-1,4) - x0po(iFam-2, 4);
%     dyp = x0po(iFam-1,5) - x0po(iFam-2, 5);
%     dzp = x0po(iFam-1,6) - x0po(iFam-2, 6);
    dx = -0.0001 ;
    x0po1 = round(x0po(iFam-1,1), 4);
    
    dt = T(iFam-1) + 0.01;
    
    % new initial guess
    x0po_g = [x0po1 + dx;  x0po(iFam-1,2); 0;...
              x0po(iFam-1,4); x0po(iFam-1,5); ...
              x0po(iFam-1,6)];
    t1po_g = (T(iFam-1) + dt)/2 + delt;
    
    % differential correction
    [x0po_iFam, t1_iFam] = poSVD3BP3d(x0po_g, mu, epsilon) ;
    
    x0po(iFam, 1:N) = x0po_iFam' ;
    T(iFam) = t1_iFam ;
    save temp_x0po x0po
    save temp_T T
end

% Store the family of initial states
dum = [x0po T];
% save x0po_T_SVD_L1_2.dat -ascii -double dum

end