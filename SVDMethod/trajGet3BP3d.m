function [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS)

global param

param = mu;

MODEL = 'cr3bp_svd';

if nargin==3
    OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12); % lower accuracy
    % OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-16); % high accuracy
end

TSPAN = linspace(0, tf, 5000);

[t,x] = ode113(MODEL, TSPAN, x0, OPTIONS);

end