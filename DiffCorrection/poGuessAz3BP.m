function Az = poGuessAz3BP(mu, eqNum, lambda, k, Ax)

%Az = poGuessAz3BP(mu, eqNum, Ax)
% 
% Generates Az based on given Ax
% This is not actually used in our code to find periodic orbits
% Usually we set Az=Ax, or Az=Ax*10, or Az=Ax/10
%
% Follow Eq (3.46)

c2 = cn3BP(mu, eqNum, 2);
c3 = cn3BP(mu, eqNum, 3);
c4 = cn3BP(mu, eqNum, 4);
d1 = 3*lambda^2*(k*(6*lambda^2 - 1) - 2*lambda)/k;
d21 = -c3/2/lambda^2;

a21 = 3*c3*(k^2 - 2)/4/(1 + 2*c2);
a22 = 3*c3/4/(1 + 2*c2);
a23 = -3*c3*lambda*(3*k^3*lambda - 6*k*(k-lambda) + 4)/(4*k*d1);
a24 = -3*c3*lambda*(3*k*lambda + 2)/(4*k*d1);

a1 = -(3/2)*c3*(2*a21 + a23 + 5*d21) - (3/8)*c4*(12 - k^2);
a2 = (3/2)*c3*(a24 - 2*a22) + (9/8)*c4;

b21 = -3*c3*lambda*(3*k*lambda - 4)/2/d1;
b22 = 3*c3*lambda/d1;

t = 1/(2*lambda*(lambda*(1+k^2) - 2*k));

s1 = ( (3/2)*c3*(2*a21*(k^2-2) - a23*(k^2+2) - 2*k*b21) -...
    (3/8)*c4*(3*k^4 - 8*k^2 + 8) )*t;
s2 = ( (3/2)*c3*(2*a22*(k^2-2) - a24*(k^2+2) - 2*k*b22 + 5*d21) +...
    (3/8)*c4*(12 - k^2) )*t;

l1 = a1 + 2*lambda^2*s1;
l2 = a2 + 2*lambda^2*s2;

Delta = lambda^2 - c2;

% Ax >= 0.1024

Az2 = -(l1*Ax^2 + Delta)/l2;

Az = sqrt(Az2);

end