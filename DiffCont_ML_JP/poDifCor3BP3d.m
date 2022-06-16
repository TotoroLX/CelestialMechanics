function [x0po, t1] = poDifCor3BP3d(x0, mu, tf)

%[x01, t1] = poDifCor3BP3d(x0poGuess1, mu)
%
% Differential correction routine to create a halo orbit
% based on Shane Ross's code
%
% inputs:
%   x0: initial guess
%   mu: mass parameter
%
% outputs:
%   x0po: initial condition after differential correction
%   t1: half-period
%
% Here we assume dt1 = 0

global param

param = mu;

% x0 = [1.034347368116448;0;-0.189464797724911;0;-0.128589215060593;0];
% set tolerances for integration and perpendicular crossing of xz-plane
MAXxdot = 1.e-10;
MAXzdot = 1.e-4;
RelTol = 2.5e-14;
AbsTol = 1.e-22;

MAXattempt = 200;

xdot1 = 1;
attempt = 0;

show = 0;
if show == 1
    figure()
    plot3(1-mu,0,0,'ro')
    hold on
end

if nargin < 3
    tf = 1.5;
end

while abs(xdot1) >= MAXxdot
    if attempt >= MAXattempt
        ERROR = 'Maximum iterations exceed';
        disp(ERROR);
        break
    end
    
    % find first xz-plane crossing
    
    MODEL = 'cr3bp';
    TSPAN = linspace(0, tf, 5000);
    % Events ON for differential correction
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events','on');
    [tt, xx, t1, xx1, i1] = ode113(MODEL, TSPAN, x0, OPTIONS);
    
    t1 = tt(end);
    x1 = xx(end,1);
    y1 = xx(end,2);
    z1 = xx(end,3);
    xdot1 = xx(end,4);
    ydot1 = xx(end,5);
    zdot1 = xx(end,6);
    
    % Compute the state transition matrix
    
    OPTIONS = odeset('RelTol', RelTol, 'AbsTol', AbsTol); % No events
    
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS);% try tf->t1
    
    attempt = attempt + 1;
    fprintf('::poDifCor : iteration %d\n',attempt);
    
    % show =1 to plot successive orbit (default = 0)
    %show = 0; % set to 1 to plot successivly
    if show == 1
        plot3(x(:,1),x(:,3),x(:,2), '.-');
        hold on
        grid on
        xlabel('X')
        ylabel('Z')
        zlabel('Y')
        plot3(x(1,1),x(1,3),x(1,2),'b*')
        plot3(x(end,1),x(end,3),x(end,2),'bo')
        pause(0.01);
    end
    
    % correction initial condition
    r1 = sqrt( (x1+mu)^2 + y1^2 + z1^2 ); % r1 distance to m1
    r2 = sqrt( (x1-(1-mu))^2 + y1^2 + z1^2 ); % r2 distance to m2

    Ux1 = x1 - (1-mu)*(x1+mu)/r1^3 - mu*(x1-(1-mu))/r2^3;
    Uz1 = -(1-mu)*z1/r1^3 - mu*z1/r2^3;

    xdotdot1 = 2*ydot1 + Ux1;  % x acceleration
    zdotdot1 = Uz1; % z acceleration
    
    %dx0 = 0
    phi_ = [phi_t1(2,3), phi_t1(2,5), ydot1; ...
        phi_t1(4,3), phi_t1(4,5), xdotdot1; ...
        phi_t1(6,3), phi_t1(6,5), zdotdot1];
    %dz0 = 0
    phi_ = [phi_t1(2,1), phi_t1(2,5), ydot1; ...
        phi_t1(4,1), phi_t1(4,5), xdotdot1; ...
        phi_t1(6,1), phi_t1(6,5), zdotdot1];
    
    b = [-y1; -xdot1; -zdot1];
    dx = inv(phi_)*b;
    
    % I added this cheat for large values of the parameter
    % based on experience with convergence problems
    if mu > 1.e-3
        DAMP = 1-0.99^attempt ;
        %DAMP = 1-0.98^attempt ;
    else
        DAMP = 1 ;
    end
    
    %tf
    %if assume dx0 = 0
    %x0(3) = x0(3) + DAMP*dx(1);
    % if assume dz0 = 0
    x0(1) = x0(1) + DAMP*dx(1);
    x0(5) = x0(5) + DAMP*dx(2);
    
    epsilon = 1;
    if abs(xdot1) < 1.e-4
        epsilon = 0;
    end
    if x(end,2)*x(1,5)>0
        tf = t1 + dx(3) + epsilon*0.001;
    elseif x(end,2)*x(1,5)<0
        tf = t1 + dx(3) - epsilon*0.001;
    end
    %tf = tf + 10; % this is only for Ax=0.0001,Az=0.001 (and chaos shows for
                  %  this initial condition)
                  
    if attempt == 6 % change to 50 for halo orbits 1
        xdot1
        tf-t1;
    end
    
end
xdot1
x0po = x0;

end