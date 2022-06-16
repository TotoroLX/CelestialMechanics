function [x0po, t1] = poSVD3BP3d(x0, mu, epsilon, tf)

%[x0po, t1] = poSVD3BP3d(x0, mu)
%
% Periodic orbits generation method by Ryan P. Russell
%
% inputs:
%   x0: initial guess
%   mu: mass parameter
%
% outputs:
%   x0po: initial condition after differential correction
%   t1: half-period

global param

param = mu;

% C = jacobiConst(x0, mu)
% d = missDist(x0, x1, C0, C1)
d = 1;

RelTol = 1.e-6;
AbsTol = 1.e-12;
MAXd = 1.e-6;

MAXattempt = 1000;

attempt = 0;

show  = 1;
if show == 1
    figure()
    hold on
    plot3(1-mu, 0, 0, 'rs')
%     plot3(-1.005062644785412, 0, 0, 'rs')
    hold on
end

if nargin < 4
    tf = 10;
end

while d > MAXd
    if attempt >= MAXattempt
        ERROR = 'Maximum iterations exceed';
        disp(ERROR);
        break
    end
    
    MODEL = 'cr3bp_svd';
    TSPAN = linspace(0, tf, 2000);
    OPTIONS = odeset('RelTol', RelTol, 'AbsTol',AbsTol,'Events', 'on');
    [tt, xx, t1, xx1, i1] = ode113(MODEL, TSPAN, x0, OPTIONS);
    
    t1 = tt(end);
    x1 = xx(end,1);
    y1 = xx(end,2);
    z1 = xx(end,3);
    xdot1 = xx(end,4);
    ydot1 = xx(end,5);
    zdot1 = xx(end,6);
    
    % compute the state transition matrix
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);
    [x,t,phi_t1, PHI] = stateTransMat3BP3d(x0, t1, mu, OPTIONS);
    
    dcdx = derivJacobi(x(end,:), mu);
    
    attempt = attempt + 1;
    
    if mod(attempt, 100) == 0
        fprintf('::poSVD : iteration %d\n', attempt);
    end
    
    if show == 1 && mod(attempt,1) == 0
        plot3(x(:,1), x(:,2), x(:,3), '.-', 'MarkerSize',2);
        hold on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        plot3(x(1,1), x(1,2), x(1,3), 'b*', 'MarkerSize',2)
        plot3(x(end,1), x(end,2), x(end,3), 'bo', 'MarkerSize',2)
        pause(0.01);
    end
    
    % singular value decomposition correction
    idx = [1 2 4 5 6];
    C0 = jacobiConst(x0, mu);
    C1 = jacobiConst(x(end,:),mu);
    dKdxi0Inv = svdCorrection(phi_t1, x(end,:), mu, dcdx, epsilon);
    K = [x(end, idx) - x0(idx)' (C1 - C0)];
    
    dx0 = -dKdxi0Inv*K(1:5)';
    
    if mu > 1.e-3
        DAMP = 1;
        %DAMP = 1; %(for L1 and L2)
%         DAMP = 1;
    end
    
    x0(idx) = x0(idx) + DAMP*dx0;
    
    if mod(attempt,1) == 0
        % 
        d;
        K(1:5)
    end
    
    %tf = t1 + 0.01;
%     if x(end,3)*x(1,6)>0
%         tf = t1 - 0.001;
%     elseif x(end,3)*x(1,6)<0
%         tf = t1 + 0.001;
%     end
    
end

x0po = x0;

end