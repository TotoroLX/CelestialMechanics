function [x0po, T] = poGet(mu, nFam, x0poT)

%[x0po, T] = poGet(mu, eqNum, x0poT)
% 
% Find periodic orbits based on given periodic orbits
%
% inputs:
%   mu: mass parameter
%   nFam: 
%   x0poT: given initial conditions of periodic orbits
%
% outputs:
%   x0po: new family of initial states of periodic orbits
%   T: correspond periods

% delt = guessed change in period between successive orbits in family
delt = -1.e-2 ; % <==== may need to be changed

N = 6; % dimension of phase space

len= size(x0poT, 1) ;

x0po = zeros(nFam+1, N);
T = zeros(nFam+1, 1);

RelTol = 2.5e-14;
AbsTol = 1.e-22;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%% find new orbits
ix0po = x0poT(1, 1:6) ;
itf = x0poT(1, end) + 0.001;

x0po(1,:) = x0poT(1:6) ;
T(1) = x0poT(end) ;

iFam = 1 ;
i = 1 ;
while i>=1%:2:len
    
    %iFam = iFam ;
    fprintf('::poFamGet : number %d\n',iFam) ;
    
    [x01, t1] = poLP3BP3d(ix0po, mu, itf) ;
    
    x0po(iFam+1, 1:N) = x01(:)';
    T(iFam+1) = t1;

    iFam = iFam + 1;
    if iFam > nFam
        break
    end
    if mod(i,10) == 0
        close all
    end
    
    ix0po = x01 ;
    itf = t1 + 0.0001 ;
    i = i + 1 ;
end

%% improve the given results
% iFam = 1 ;
% for i=1:2:len
%     ix0po = x0poT(i, 1:6) ;
%     itf = x0poT(i, end) ;
%     %iFam = iFam ;
%     fprintf('::poFamGet : number %d\n',iFam) ;
%     
%     [x01, t1] = poLP3BP3d(ix0po, mu, itf) ;
%     
%     x0po(iFam, 1:N) = x01(:)';
%     T(iFam) = t1;
% 
%     iFam = iFam + 1;
%     if iFam > nFam
%         break
%     end
%     if mod(i,20) == 0
%         close all
%     end
% end

% Store the family of initial states
dum = [x0po T];
% save x0po_T.dat -ascii -double dum

end
