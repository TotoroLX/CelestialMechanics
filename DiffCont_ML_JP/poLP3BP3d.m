function [x0po, t1] = poLP3BP3d(x0, mu, tf)

%[x01, t1] = poLP3BP3d(x0poGuess1, mu)
%
% Algorithm in paper LP
%
% inputs:
%   x0: initial guess
%   mu: mass parameter
%
% outputs:
%   x0po: initial condition after differential correction
%   t1: half-period
% 
% Models: CR3BP
% dot{X} = F(X(t),t)
%
% Here we assume dt1 = 0

% set tolerances for integration
RelTol = 2.5e-14;
AbsTol = 1.e-22;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);


MAXattempt = 1500;
attempt = 1;

% variation of Jacobi constant
dc = 0.1 ;

show = 1;
if show == 1
    figure()
    hold on
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    plot3(x(:,1),x(:,2),x(:,3),'g.', 'MarkerSize', 2);
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    %plot3(0.836915146106164,0,0,'rs')
    plot3(1-mu,0,0,'ro')
    %plot3(-mu,0,0,'ro')
end

% initially, the initial states are equal to final states
% note usually we use xi as the initial condition
xi = x0;

deltaxiNorm = 1e10;
maxdeltaxi = 1e10;

Tol = 1e-11;

errors = [];

while attempt <= MAXattempt
    
    
    % now we also need stateTransMat due to eqns (6,7)
    [x,t,phi_tf,PHI] = stateTransMat3BP3d(xi, tf, mu, OPTIONS);
    phi = phi_tf - diag([1,1,1,1,1,1]);
    xT = x(end, :) ;
    
    if show == 1 && mod(attempt, 20) == 0
        plot3(x(:,1),x(:,2),x(:,3), '.-', 'MarkerSize', 2);
        hold on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2)
        plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2)
        pause(0.01);
    end
    
    % see equation (6,7) F
    dotX = cr3bpF(mu, xi) ;
    
    % equation * in my notes
    dCdxi = derivJacobi(xi, mu);
    
    mat = [dotX, phi; 0, dCdxi];
    
    if attempt >= 2
        dc = 0;
    end    
    
    deltaxi = xi - xT ;
    dtdxi = inv(mat)*[deltaxi';dc] ;
    if max(abs(dtdxi)) > 0.001
        DAMP = 0.01;
        dtdxi = DAMP.* dtdxi ;
%         if attempt > 750
%             DAMP = 0.001;
%             dtdxi = DAMP.* dtdxi ;
%         end
    else
        DAMP = 1;
    end
    dt = dtdxi(1) ;
    dxi = dtdxi(2:7) ;
    
    errorN = norm(deltaxi(1:3))/norm(xi(1:3)) +...
            norm(deltaxi(4:6))/norm(xi(4:6));
    errors(attempt) = errorN ;
    % record the best result
    if errorN <= deltaxiNorm
        fprintf('::poDifCor : iteration %d, norm : %.15f\n',attempt, deltaxiNorm);
        deltaxiNorm = errorN ;
        
        maxdeltaxi = max(abs(deltaxi));
        x0po = xi ;
        t1 = tf ;
        if deltaxiNorm <= Tol
            break
        end
    end
    
    % display some information
    if mod(attempt, 20) == 0
        fprintf('::poDifCor : iteration %d, norm : %.15f\n',attempt, errorN);
    end
    
    
    attempt = attempt + 1;
    % update the initial condition and the period
    xi = xi + dxi';
    tf = tf + dt;
end

ERROR = 'Maximum iterations exceed';
disp(ERROR);

figure();plot(errors)

end