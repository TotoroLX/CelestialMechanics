% some test
muM = 4902.799; % km3/s2 -- mu Moon from Vallado
muE = 3.986004415e5; % km3/s2 -- mu Earth from Vallado
mu = muM/(muE + muM); % 0.012150581477177

% mu=3.054248395728148e-06;

eqNum = 1;

%x0poT = x0po_T_allDCL3 ;
%x0poT = sortrows(x0poT, 1, 'descend');
x0poT = [ -1.000722402232119  -0.122704917861421                   0  -0.223378755683497   1.812129142010293   0.566515929039525  6.311079455031661 ];

nFam = 5; %size(x0poT, 1) ;

% Used differential correction
[x0po, T] = poGet(mu, nFam, x0poT);

%RelTol = 3.e-14 ; AbsTol = 1.e-16; % high accuracy;
%RelTol = 3.e-10 ; AbsTol = 1.e-12; % low accuracy
%RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy
RelTol = 2.5e-14;
AbsTol = 1.e-22;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);


figure()
hold on
grid on

ekk = [];
err = zeros( size(x0po, 1), 1);

for k = 1:size(x0po,1)-14%nFam
    x0 = x0po(k,:) ;
    tf = T(k) ;
    %[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    [x,t,phi_tf,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS);
    eigs = eig(phi_tf) ;
    ekk(k,:) = eigs ;
    err(k) = max( abs( x(end,1:3) - x0(1:3) ) ) ;
    plot3(x(:,1),x(:,2),x(:,3),'.-', 'MarkerSize', 2, 'MarkerEdgeColor','#D95319');
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    hold on
    pause(0.0001) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3(0.836915146106164,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(1.155733835140644,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(-1.005062644785412,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(1-mu,0,0,'ko', 'MarkerFaceColor', 'r')
plot3(-mu,0,0,'ko', 'MarkerFaceColor', 'r')
