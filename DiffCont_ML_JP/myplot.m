function myplot(x0,tf,mu)

RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy

OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

figure()
hold on
[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
plot3(x(:,1),x(:,3),x(:,2),'b.-');
plot3(x(1,1),x(1,3),x(1,2),'r*');
plot3(x(end,1),x(end,3),x(end,2),'ro');

xlabel('X');
ylabel('Z');
zlabel('Y');
grid on
plot3(0.836915146106164,0,0,'rs')
plot3(1-mu,0,0,'ro')

    
end