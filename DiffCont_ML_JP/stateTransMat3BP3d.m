function [x,t,phi_tf,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS)

%[x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, t1, mu, OPTIONS)
%
% Computes the state transition matrix, phi_tf = PHI(0,tf)
% based on Shane Ross
%
% inputs:
%   x0: initial state
%   tf: half-period of from initial integration
%   mu: mass parameter
%
% outputs:
%   x: states from t=0 to t=tf
%   t: time corresponding to states x
%   phi_t1: phi(t1)
%   PHI: state transition matrix phi

global param

param = mu;

N = length(x0);

MODEL = 'sysSolveCRTBP';   % variational equation for 3BP
MODEL = 'vareqn3bp';   % variational equation for 3BP

TSPAN = linspace(0, tf, 15000); % time interval based on guess integration

PHI_0(1:N^2) = reshape(eye(N),N^2,1); % initial condition
PHI_0(1+N^2:N+N^2) = x0; % phi should be simultaneously integrated with states

[t, PHI] = ode113(MODEL, TSPAN, PHI_0, OPTIONS,[],mu);
%[t, PHI] = ode113(MODEL, [0 tf], PHI_0, OPTIONS,[],mu);

x = PHI(:, 1+N^2:N+N^2); % trajectory from 0 to tf

c = PHI(end,1:N^2);

phi_tf = reshape(PHI(length(t), 1:N^2), N, N); %phi(0,tf)

for i=1:N
    for j=1:N
        phi_tf(i,j) = c(N*(i-1)+j);
    end
end

end
