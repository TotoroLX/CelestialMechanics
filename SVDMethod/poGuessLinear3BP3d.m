function [x0poGuess, TGuess] = poGuessLinear3BP3d(mu, eqNum, Ax)

%[x0poGuess, TGuess] = poGuessLinear3BP(mu, eqNum, Ax)
%
% Uses linearised model around the equilibrium point to get a first guess
% for generating a periodic orbit, assume z = 0
%
% inputs:
%   mu: mass parameter
%   eqNum: the number of eq point
%   Ax: x-amplitude of periodic orbit
%
% outputs:
%   x0poGuess: initial guessed state for halo orbit
%              [x0 y0 0 0 0 yv0]
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

% Az = poGuessAz3BP(mu, eqNum, lambda, k, Ax); % not used
Ay = Ax*10;%
% Az = abs(Az);
%Az = 0;

x0poGuess = zeros(6,1);
x0poGuess(1) = eq(1) - Ax;
x0poGuess(2) = eq(2) + Ay*sin(1*pi/2);   % we first test n=1
x0poGuess(6) = eq(6) + Ax*k*lambda;

% use the point we get before
% and slightly modify it
x0poGuess = [0.890;0.001;0;...
    0;0;0.01];
x0poGuess = [.968531905820430;...
    -0.016539104132630;...
    0;...
    0.531354597072870;...
    -0.568790585213285;...
    0.167349497962551];
x0poGuess = [0.988360555400675;0.004481176568162;0;...
    -0.098106749679581;-1.644397521594456;1.640416175720395];
if Ax == 0.02
    x0poGuess = [0.844765770471713;-0.056805336483208;0;...
        -0.035952982443933;-0.014098375874756;0.025389558388583];
else
    x0poGuess = [0.844647937718099;-0.056447317105796;0;...
        -0.035690452273597;-0.014094362504229;0.021170335874687];
end

TGuess = 2*pi/lambda;

end