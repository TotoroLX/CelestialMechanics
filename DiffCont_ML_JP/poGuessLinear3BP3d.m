function [x0poGuess, TGuess] = poGuessLinear3BP3d(mu, eqNum, Ax)

%[x0poGuess, TGuess] = poGuessLinear3BP(mu, eqNum, Ax)
%
% Uses linearised model around the equilibrium point to get a first guess
% for generating a periodic orbit
%
% inputs:
%   mu: mass parameter
%   eqNum: the number of eq point
%   Ax: x-amplitude of periodic orbit
%
% outputs:
%   x0poGuess: initial guessed state for halo orbit
%              [x0 0 z0 0 yv0 0]
%   TGuess: periodic of periodic orbit

eq = eqPoint3BP(mu, eqNum);   % location of eq point
eq = [eq 0 0 0]; % L1:0.836915146106164 L2:1.155733835140644
                 % L3:-1.005062644785412
% Get the Jacobian matrix and its eigenvalues
DU = DUmatrix3BP(eq, mu);
eqPointEig = eig(DU);
[M,~] = size(DU);

for i=1:M-2
    if abs(real(eqPointEig(i))) < 1e-10
        lambda = abs(imag(eqPointEig(i))); % lambda
        break
    end
end

c2 = cn3BP(mu, eqNum, 2);

k = (lambda^2 + 1 + 2*c2)/(2*lambda);

Az = poGuessAz3BP(mu, eqNum, lambda, k, Ax); % not used
Az = Ax/10;%
% Az = abs(Az);
%Az = 0;

x0poGuess = zeros(6,1);
x0poGuess(1) = eq(1) - Ax;
x0poGuess(3) = eq(3) + Az*sin(1*pi/2);   % we first test n=1
x0poGuess(5) = eq(5) + Ax*k*lambda;

% x0poGuess(1) = 0.836815146106164;
% x0poGuess(3) = -1.279390117323268;
% x0poGuess(5) = -0.616836892970480;

TGuess = 2*pi/lambda;

end