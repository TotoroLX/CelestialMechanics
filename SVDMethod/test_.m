% some test
muM = 4902.799; % km3/s2 -- mu Moon from Vallado
muE = 3.986004415e5; % km3/s2 -- mu Earth from Vallado
mu = muM/(muE + muM); % 0.012150581477177

epsilon = 1.e-4;


%RelTol = 3.e-14 ; AbsTol = 1.e-16; % high accuracy;
%RelTol = 3.e-10 ; AbsTol = 1.e-12; % low accuracy
RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy

OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%%

x0_=[[0.943553185939109;0.015497025068962;0;...
    0.531332017860071;0.123268245892669;0.492528079554908;2.626258844894558],...
    [0.943044654344087;0.015457998465568;0;...
    0.529230173375609;0.126443076098463;0.489184926164547;2.631633602321876],...
    [0.982856786261038;0.006870087813461;0;...
    0.811403580811102;-0.825753551531886;1.235426483336298;2.162928613406002],...
    [0.981663286978653;0.007315903443993;0;...
    0.835888084370029;-0.699097787821255;1.163991740515774;2.178107511735337],...
    [1.005479389865741;0.013646216133220;0;...
    -0.495538562327786;-0.639498716376086;0.665866054417634;2.263995684068486],...
    [1.005052153879725;0.012903125990372;0;...
    -0.517336624975045;-0.638267215803521;0.681107898103560;2.266515810647406],...
    [1.005831296108898;0.013219448029302;0;...
    -0.510967053044670;-0.626205244003947;0.666891491237200;2.273570879725431],...
    [1.007397738688967;0.013878844597625;0;...
    -0.497932060595133;-0.604614567048500;0.640203339837680;2.288103871185747],...
    [1.016359937117348;0.017471419518325;0;...
    -0.434112286607389;-0.521287544848413;0.528583777389106;2.368423122606387],...
    [1.031221934685759;0.024010801655346;0;...
    -0.352089065107875;-0.457977577671954;0.418935108976313;2.501564132191694],...
    [1.026953597092422;0.022361440096865;0;...
    -0.369253965977198;-0.472923197320486;0.443083835108287;2.461519765759422],...
    [1.030886158712089;0.024166166515731;0;...
    -0.350314110367999;-0.460506443289485;0.419640575009179;2.496125236219684],...
    [1.009271330859460;0.015333412320826;0;...
    -0.466288176931208;-0.590286317286946;0.604999021698881;2.299382820723760],...
    [1.057726010973495;0.038554376818750;0;...
    -0.250681335420820;-0.419185343065547;0.314751047117090;2.741597991432257],...
    [1.064738587680264;0.043271849989847;0;...
    -0.228972884865106;-0.415519399604692;0.296458227967495;2.806484236626749],...
    [1.141757375396668;0.161461289991517;0;...
    -0.003266212653614;-0.435288086168592;0.184094688871538;3.702827564582187],...
    [1.1119594054551374;0.0;0;...
    0.0;-0.18124841270866476;0.43576411644659757; 4.423308886403689],...
    [0.87642891534439826;0.0;0.0;...
    0.0;0.0090290217274574605;0.46553952107519370; 4.060839820175107]];

datFiles1 = dir("x0po_T_SVD_L1*.dat");
N1 = length(datFiles1);% 1020 NRO
x0po1 = [];
for i=1:N1
    thisFile = datFiles1(i).name ;
    x0poi = load(thisFile) ;
    x0po1 = [x0po1; x0poi] ;
end


datFiles2 = dir("x0po_T_SVD_L2*.dat");
N2 = length(datFiles2);
x0po2 = [];
for i=1:N2
    thisFile = datFiles2(i).name ;
    x0poi = load(thisFile) ;
    x0po2 = [x0po2; x0poi] ;
end


datFiles3 = dir("x0po_T_SVD_L3*.dat");
N3 = length(datFiles3);
x0po3 = [];
for i=1:N3
    thisFile = datFiles3(i).name ;
    x0poi = load(thisFile) ;
    x0po3 = [x0po3; x0poi] ;
end

x0po = [x0po1; x0po2; x0po3] ;

x0po_iter = x0po1;

x0po_iter = sortrows(x0po_iter, 1);

% % select orbits
idxx = [1020] ;
%difference = [];
j=2;
for k = 1:size(x0po_iter,1)
    x0 = x0po_iter(k, 1:6) ;
    tf = x0po_iter(k, end);
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    %difference(k)= max(abs(x(end,:) - x0));
    if max(abs(x(end,:) - x0)) < 1.e-6
        idxx(j) = k ;
        j = j+1 ;
    end
end
x0po_iter = x0po_iter(idxx,:);
% %

x0po_iter = clearx0(x0po_iter, mu, OPTIONS);

diffs = zeros(size(x0po_iter, 1), 1) ;

figure()
hold on
for k = 1:9:size(x0po_iter, 1)
    x0 = x0po_iter(k, 1:6) ;
    tf = x0po_iter(k, end);
    %[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    diffs(k) = max( abs( x0(1:3) - x(end,1:3) ) ) ;
    plot3(x(:,1),x(:,2),x(:,3),'b.-', 'MarkerSize',2);
    hold on
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize',2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize',2);
    
    pause(0.001) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3(0.836915146106164,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(1.155733835140644,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(-1.005062644785412,0,0,'ks', 'MarkerFaceColor', 'r')
plot3(1-mu,0,0,'ko', 'MarkerFaceColor', 'r')
plot3(-mu,0,0,'ko', 'MarkerFaceColor', 'r')
grid on

title("Periodic orbits with SVD method")


% % %
x0 = [0.9292, 0, 0.2914, 0, 0.0817, 0] ;
tf = 2.1509 ;
%x0 = [1.0118, 0, 0.1739, 0, -0.0799, 0] ;
%tf = 1.3743 ;
%x0 = [0.987615013885612, 0, -0.004758777893685, 0, 2.234704092755416, 0] ;
%tf = 0.687098687215876*2 ;
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
err = zeros(len, 1);

figure()
hold on
for k = 1:len
    x0 = x0po_iter(k,1:6) ;
    tf = x0po_iter(k, end) ;
    %[x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    ek = eig(phi_t1) ;
    ekk(k,:) = ek' ;
    err(k) = max( abs(x(end,1:3) - x0(1:3)) );
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

%% JacobiConstant and Perilune and period

len = size(x0po_iter, 1) ;
Jc = zeros(len, 1);
perilune = zeros(len, 1);

for idx = 1:len
    Jc(idx) = jacobiConst(x0po_iter(idx,:), mu);
    
    x0 = x0po_iter(idx,1:6) ;
    tf = x0po_iter(idx, end) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    x = x(:,1:3) - [1-mu,0,0];
    perilune(idx) = min( vecnorm(x, 2, 2) );
end
figure()
plot(x0po_iter(1:end, 1), Jc(1:end), 'bo', 'MarkerSize', 2)
title("Jacobi Constant vs x-coordinate")
xlabel("X")
ylabel("Jacobi Constant")
figure()
plot(x0po_iter(1:end, 2), Jc(1:end), 'bo', 'MarkerSize', 2)
title("Jacobi Constant vs y-coordinate")
xlabel("Y")
ylabel("Jacobi Constant")
figure()
plot(perilune(1:end), Jc(1:end), 'bo', 'MarkerSize', 2)
title("Jacobi Constant vs Perilune Radius")
xlabel("Perilune Radius")
ylabel("Jacobi Constant")
figure()
plot(perilune(1:end), x0po_iter(1:end, end), 'bo', 'MarkerSize', 2)
title("Period vs Perilune Radius")
xlabel("Perilune Radius")
ylabel("Period")
figure()
plot(Jc(1:end), x0po_iter(1:end, end), 'bo', 'MarkerSize', 2)
title("Period vs Jacobi Constant")
ylabel("Period")
xlabel("Jacobi Constant")

%% ignore this
temp_x0po_iter = zeros(100, 7) ;

level0 = 1.1137;
idx = 1;

for k = 1:size(x0po_iter,1)
    x0 = x0po_iter(k, 1:6) ;
    if x0(1) > level0
        temp_x0po_iter(idx, :) = x0po_iter(k,:) ;
        level0 = level0 + 0.0022 ;
        idx = idx + 1 ;
    end
end