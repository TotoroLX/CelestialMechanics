function cn = cn3BP(mu, eqNum, n)

%cn = cn3BP(mu, eqNum, n)
%
% Computes some parameters needed for computation
% see Appendix A
%
% inputs:
%   mu: mass parameter
%   eqNum: the number of eq point
%   n: index of cn parameter
%
% outputs:
%   cn: parameter of interest
%
% See Appendix A, [25] Eq (43 -45)

alpha = (mu/3/(1-mu))^(1/3);

gamma1 = alpha*(1 - alpha/3 - alpha^2/9 - 23*alpha^3/81);
gamma1 = gamma3BP(mu,1);
gamma2 = alpha*(1 + alpha/3 - alpha^2/9 - 31*alpha^3/81);
gamma3 = 1 - 7*mu/12*(1 + (7*mu/12)^2*23/84);

cn=0; % for l4 and l5

if eqNum == 1
    cn = (mu + (-1)^n*(1-mu)*gamma1^(n+1)/(1-gamma1)^(n+1))/gamma1^3;
elseif eqNum == 2
    cn = (-1)^n*(mu + (1-mu)*gamma2^(n+1)/(1+gamma2)^(n+1))/gamma2^3;
elseif eqNum == 3
    cn = (1-mu + mu*(gamma3/(1+gamma3))^(n+1))/gamma3^3;
end

end