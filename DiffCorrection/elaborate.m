% This file will try to increase the accuracy of initial conditions
%
% Assume the velocity is perpendicular to the xy-plane
%
muM = 4902.799; % km3/s2 -- mu Moon from Vallado
muE = 3.986004415e5; % km3/s2 -- mu Earth from Vallado
mu = muM/(muE + muM); % 0.012150581477177

RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy

OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);


RelTol = 2.5e-14;
AbsTol = 1.e-22;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%% load all the saved data and run the code
datFiles1 = dir("x0po_T_L1*.dat");
N1 = length(datFiles1);
x0po1 = [];
for i=1:N1
    thisFile = datFiles1(i).name ;
    x0poi = load(thisFile) ;
    x0po1 = [x0po1; x0poi] ;
end


datFiles2 = dir("x0po_T_L2*.dat");
N2 = length(datFiles2);

x0po2 = [];
for i=1:N2
    thisFile = datFiles2(i).name ;
    x0poi = load(thisFile) ;
    x0po2 = [x0po2; x0poi] ;
end


datFiles3 = dir("x0po_T_L3*.dat");
N3 = length(datFiles3);
x0po3 = [];
for i=1:N3
    thisFile = datFiles3(i).name ;
    x0poi = load(thisFile) ;
    x0po3 = [x0po3; x0poi] ;
end

x0po = [x0po1; x0po2; x0po3] ;

x0po_iter = x0po1;

x0po_iter = sortrows(x0po_iter,1);

len = size(x0po_iter, 1) ;
x0po = zeros(len, 7) ;
for i = 1:len
    x0 = x0po_iter(i, 1:6) ;
    tf = x0po_iter(i, 7) ;
    fprintf('::poFam : number %d\n',i) ;
    [nx0, ntf] = npoDifCor3BP3d(x0, mu) ;
    x0po(i,:) = [nx0, ntf] ;
end

x0po_iter = clearx0(x0po_iter, mu, OPTIONS);

figure()
hold on
grid on
for k = 1:size(x0po_iter,1)
    x0 = x0po_iter(k,1:6) ;
    tf = x0po_iter(k, end) ;
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    max( abs(x(end,:)-x0 ))
    plot3(x(:,1),x(:,2),x(:,3),'b.-', 'MarkerSize', 2);
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    hold on
    pause(0.01) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3(0.836915146106164,0,0,'ks')
plot3(1.155733835140644,0,0,'ks')
plot3(-1.005062644785412,0,0,'ks')
plot3(1-mu,0,0,'ko')
plot3(-mu,0,0,'ko')
title("Halo orbits, Differential Correction")

% % %
%x0 = [0.9292, 0, 0.2914, 0, 0.0817, 0] ;
%tf = 2.1509 ;
x0 = [1.0118, 0, 0.1739, 0, -0.0799, 0] ;
tf = 1.3743 ;
x0 = [0.987615013885612, 0, -0.004758777893685, 0, 2.234704092755416, 0] ;
tf = 0.687098687215876*2 ;
[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
plot3(x(:,1),x(:,2),x(:,3),'g.-', 'MarkerSize', 2);
plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
% % %
%% monodromy matrix

len = size(x0po_iter, 1) ;
lambda1 = zeros(len,2) ;
lambda2 = zeros(len,2) ;
lambda3 = zeros(len,2) ;
lambda4 = zeros(len,2) ;
lambda5 = zeros(len,2) ;
lambda6 = zeros(len,2) ;
ekk = zeros(len, 6);

figure()
hold on
for k = 1:len
    x0 = x0po_iter(k,1:6) ;
    tf = x0po_iter(k, end) ;
    %[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    ek = eig(phi_t1) ;
    ekk(k,:) = ek' ;
    plot(ek, 'o', 'MarkerSize', 2) ;
    hold on
end
title("Eigenvalues of Monodromy matrix")
figure()
plot(ekk(:,1:2), 'o', 'MarkerSize', 2) ;
figure()
plot(ekk(:,3:4), 'o', 'MarkerSize', 2) ;
figure()
plot(ekk(:,5:6), 'o', 'MarkerSize', 2) ;
%print(gcf, 'eig11.png', '-dpng', '-r600');

%% Bifurcation
absekk = abs(ekk) ;
x = x0po_iter(:,1);
figure()
hold on
plot(x, absekk(:,3), '*', 'MarkerSize', 3)
plot(x, absekk(:,4), 's', 'MarkerSize', 3)
plot(x, absekk(:,5), 'v', 'MarkerSize', 3)
plot(x, absekk(:,6), 'o', 'MarkerSize', 3)
title("Bifurcation map")
xlabel("x")
ylabel("Eigenvalues")

