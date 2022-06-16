function dotX = cr3bpF(mu, x)
  
    m1 = 1-mu ; % Larger mass at -mu
    m2 = mu ;   % smaller mass at 1-mu
  
    r1 = sqrt( (x(1) + m2)^2 + x(2)^2 + x(3)^2 ) ;
    r2 = sqrt( (x(1) - m1)^2 + x(2)^2 + x(3)^2 ) ;
  
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
    
    dotX = xdot;
end
