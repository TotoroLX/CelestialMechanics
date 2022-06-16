function [x0po, t1] = poDifCor3BP3d(x0, mu)

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

global param

param = mu;

x0 = [0.837112192529942;0;0.001;0; 0.015537673947833;0];
x0 = [0.991975555377273;0;-0.001871657754011;0; -0.011750613711503;0];
% set tolerances for integration and perpendicular crossing of xz-plane
MAXxdot = 1.e-4;
MAXzdot = 1.e-4;
RelTol = 3.e-12;
AbsTol = 1.e-15;

MAXattempt = 25;

xdot1 = 1;
attempt = 0;

show = 1;
if show == 1
    figure()
    %plot(
end

tf = 1000;

while abs(xdot1) > MAXxdot
    if attempt > MAXattempt
        ERROR = 'Maximum iterations exceed';
        disp(ERROR);
        break
    end
    
    % find first xz-plane crossing
    
    MODEL = 'cr3bp';
    TSPAN = [0 tf];
    % Events ON for differential correction
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events','on');
    [tt, xx, t1, xx1, i1] = ode113(MODEL, TSPAN, x0, OPTIONS);
    
    if ~isempty(t1)
        t1=t1(end);
        xx1 = xx1(end,:);
        i1 = i1(end);
        x1 = xx1(1);
        y1 = xx1(2);
        z1 = xx1(3);
        xdot1 = xx1(4);
        ydot1 = xx1(5);
        zdot1 = xx1(6);
    else
        t1 = tt(end);
        x1 = xx(end,1);
        y1 = xx(end,2);
        z1 = xx(end,3);
        xdot1 = xx(end,4);
        ydot1 = xx(end,5);
        zdot1 = xx(end,6);
    end
    
    % Compute the state transition matrix
    
    OPTIONS = odeset('RelTol', RelTol, 'AbsTol', AbsTol); % No events
    
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, t1, mu, OPTIONS);
    
    attempt = attempt + 1;
    fprintf('::poDifCor : iteration %d\n',attempt);
    
    % show =1 to plot successive orbit (default = 0)
    %show = 1; % set to 1 to plot successivly
    if show == 1
        plot3(x(:,1),x(:,2),x(:,3), '.-');
        hold on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        plot3(x(1,1),x(1,2),x(1,3),'b*')
        plot3(x(end,1),x(end,2),x(end,3),'bo')
        pause(0.01);
    end
    
    % correction initial condition
    r1 = sqrt( (x1+mu)^2 + y1^2 + z1^2 ); % r1 distance to m1
    r2 = sqrt( (x1-(1-mu))^2 + y1^2 + z1^2 ); % r2 distance to m2

    Ux1 = x1 - (1-mu)*(x1+mu)/r1^3 - mu*(x1-(1-mu))/r2^3;
    Uy1 = y1 - (1-mu)*y1/r1^3 - mu*y1/r2^3;
    Uz1 = -(1-mu)*z1/r1^3 - mu*z1/r2^3;

    xdotdot1 = 2*ydot1 + Ux1;  % x acceleration
    zdotdot1 = Uz1; % z acceleration
    
    Axdot1 = phi_t1(4,1) - phi_t1(2,1)*xdotdot1/ydot1;
    Bxdot1 = phi_t1(4,3) - phi_t1(2,3)*xdotdot1/ydot1;
    Cxdot1 = phi_t1(4,5) - phi_t1(2,5)*xdotdot1/ydot1;
    Azdot1 = phi_t1(6,1) - phi_t1(2,1)*zdotdot1/ydot1;
    Bzdot1 = phi_t1(6,3) - phi_t1(2,3)*zdotdot1/ydot1;
    Czdot1 = phi_t1(6,5) - phi_t1(2,5)*zdotdot1/ydot1;
    
    % if assume dx0 = 0
    %dz0 = (zdot1*Cxdot1 - xdot1*Czdot1)/(Bxdot1*Czdot1 - Bzdot1*Cxdot1);
    %dydot0 = (xdot1*Bzdot1 - zdot1*Bxdot1)/(Bxdot1*Czdot1 - Bzdot1*Cxdot1);
    %dt1 = -(phi_t1(2,3)*dz0 + phi_t1(2,5)*dydot0)/ydot1
    
    % if assume dz0 = 0
    dx0 = (zdot1*Cxdot1 - xdot1*Czdot1)/(Axdot1*Czdot1 - Azdot1*Cxdot1);
    dydot0 = (xdot1*Azdot1 - zdot1*Axdot1)/(Bxdot1*Czdot1 - Bzdot1*Cxdot1);
    dt1 = -(phi_t1(2,1)*dx0 + phi_t1(2,5)*dydot0)/ydot1
    
    % I added this cheat for large values of the parameter
    % based on experience with convergence problems
    if mu > 1.e-3
        DAMP = 1-0.5^attempt ;
        %DAMP = 1-0.98^attempt ;
    else
        DAMP = 1 ;
    end
    
    x0
    %if assume dx0 = 0
    %x0(3) = x0(3) + DAMP*dz0;
    % if assume dz0 = 0
    x0(1) = x0(1) + DAMP*dx0;
    x0(5) = x0(5) + DAMP*dydot0;
    
    %tf = t1 + dt1;
    
end

x0po = x0;

end