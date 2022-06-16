function [out1, out2, out3] = cr3bp_svd(t, x, flag)

%[out1, out2, out3] = cr3bp(t, x, flag)
%
% model of circular restricted 3-body problem

global param

mu = param;

if nargin < 3 | isempty(flag)
    m1 = 1-mu; % LARGER MASS, at -mu
    m2 = mu; % smaller mass, at 1-mu
    
    r1 = sqrt( (x(1)+m2)^2 + x(2)^2 + x(3)^2 ); % distance to m1
    r2 = sqrt( (x(1)-m1)^2 + x(2)^2 + x(3)^2 ); % distance to m2
    
    Ux = x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-(1-mu))/r2^3;
    Uy = x(2) - (1-mu)*x(2)/r1^3 - mu*x(2)/r2^3;
    Uz = -(1-mu)*x(3)/r1^3 - mu*x(3)/r2^3;
    
    xdot = zeros(6,1);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = 2*x(5) + Ux;
    xdot(5) = -2*x(4) + Uy;
    xdot(6) = Uz;
    
    out1 = xdot;
else
    switch (flag)
        case 'events'
            if abs(t) > 0.3 % wait a short time before checking
                isterminal = 1;
            else
                isterminal = 0;
            end
            
            direction = 0; % for NRO 0, otherwise 1
            
            % check for xy-plane crossing (x(3) = 0)
            
            out1 = x(3);
            out2 = isterminal;
            out3 = direction;
    end
end

end