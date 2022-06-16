function eqPos = eqPoint3BP(mu, eqNum)

% eqPos = eqPoint3BP(mu, eqNum)
%
% inputs:
%   mu: mass parameter
%   eqNum: the number of eq point
%
% outputs:
%   eqPos: position of eq point
%
%-----------------------------------------------------------
% CR3BP CONVENTION:
%         L4
%
%
% L3-----M1-------L1---M2---L2     M1=1-mu, M2=mu
%       -mu           1-mu
%
%         L5
%
% Follow Ch 5.8.7 Prof. Casotto

alpha = (mu/3/(1-mu))^(1/3);

gamma1 = alpha*(1 - alpha/3 - alpha^2/9 - 23*alpha^3/81);
gamma1 = gamma3BP(mu,1);
l1 = 1 - mu - gamma1;
L1 = [l1 0 0];

gamma2 = alpha*(1 + alpha/3 - alpha^2/9 - 31*alpha^3/81);
l2 = 1 - mu + gamma2;
L2 = [l2 0 0];

gamma3 = 1 - 7*mu/12*(1 + (7*mu/12)^2*23/84);
l3 = -mu - gamma3;
L3 = [l3 0 0];

l4x = .5 - mu;
l4y = .5*sqrt(3);
L4 = [l4x l4y 0];

l5x = .5 - mu;
l5y = .5*sqrt(3);
L5 = [l5x l5y 0];

if eqNum == 1
    eqPos = L1;
elseif eqNum == 2
    eqPos = L2;
elseif eqNum == 3
    eqPos = L3;
elseif eqNum == 4
    eqPos = L4;
elseif eqNum == 5
    eqPos = L5;
end

end