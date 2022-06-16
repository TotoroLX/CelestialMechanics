
%
% Assume the velocity is perpendicular to the xy-plane
%
muM = 4902.799; % km3/s2 -- mu Moon from Vallado
muE = 3.986004415e5; % km3/s2 -- mu Earth from Vallado
mu = muM/(muE + muM); % 0.012150581477177

% mu=3.054248395728148e-06;

eqNum = 1;
nFam = 6;

% around l1
Ax1 = -0.15;%0.150834272416660;
Ax2 = 0.01;

%------------------------------ Halo Orbits 1 -----------------------------
%------------------------------- Halo Orbits ------------------------------
%------------------------------ Halo Orbits 2 -----------------------------
%------------------------------- Halo Orbits ------------------------------
%------------------------------- Strange Orbit ----------------------------
%x0 =[0.837215496904172;0.001000000003724;0;0;0;
%     -0.025079013876637]; t1=1.383286204771840
% {x0 = [0.875238086872289;0.011000000000000;0;0;0;
%   0.718937741293183]; t1=2.936343343809323}
% x0poGuess1=[0.823386384377028;0.01;0;0;0;0.5]
%------------------------------- Quasi Orbit -----  -----------------------
% {x0 = [0.950176635780486;-0.012000000000000;0;0;0;
%   -0.557864808600005]; t1=0.685862515063730}
% {x0 = [0.983191617799235;-0.040000000000000;0;0;0;
%   -0.549240602248344]; }
% {x0 = [0.983191617786556; 0.040000000000000;0;0;0;
%   -0.549240602240788]; }
% {x0 = [0.954392378951353;0.080000000000000;0;0;0;
%   -0.370536492666634]; }
% x0poGuess1=[0.990290839334460;0.06;0;0;0;-0.5];
% x0 = [0.974223135535817;0.06;0;0;0;0.443887077063263]
% x0=[0.823389831024856;0;0.023767681168882;0;0.135209563152812;0]; new
% x0=[0.823389831024856;0;0.003724454349521;0;0.126557411688950;0];
% x0=[0.830000000000000;0;0.114062665673766;0;0.229389349189050;0]; 2*1.393664469459206
% 0.800 to 0.844 and 0.830 to 0.933
% x0=[0.830788918804415;0;0.118926583580229;0;0.233887687281693;0]; 2*1.393102370145217;
% x0=[0.840000000000000;0;0.157513544811492;0;0.261303731572755;0]; 2*1.353475330082072;
% x0=[0.850000000000000;0;0.175465009915553;0;0.262898016243887;0]; 2*1.277199731148203;
% x0=[0.860000000000000;0;0.184070010537236;0;0.253957947138755;0]; 2*1.199815246122309;
% x0=[0.870000000000000;0;0.189332483982820;0;0.240311903325813;0]; 2*1.128235402654377;
% x0=[0.880000000000000;0;0.193280382628371;0;0.223747649523999;0]; 2*1.063769473130137;
% x0=[0.890000000000000;0;0.196839664636369;0;0.204924283340109;0]; 2*1.007208034235878;
% x0=[0.910000000000000;0;0.205477542413230;0;0.161193515502119;0]; 2*0.923458336937656;
% x0=[0.920000000000000;0;0.212788140451770;0;0.135923838975717;0]; 2*0.903153755887798;
% x0=[0.930000000000000;0;0.228428474390061;0;0.106206331514581;0]; 2*0.916717279670384;
% x0=[0.810000000000000;0;0.559469508798996;0;0.185767772787611;0]; 2*1.425794999443345;
% x0=[0.800000000000000;0;0.574135940743903;0;0.195502355356865;0]; 2*1.434404210311786;
% x0=[0.780000000000000;0;0.601824214657780;0;0.215046614649965;0]; 2*1.449157422386820;
% x0=[0.770000000000000;0;0.614927244911848;0;0.224849083100288;0]; 2*1.455524398007910;
% x0=[0.990000000000000;0;0.038540932119960;0;0.750779123869441;0]; 2*1.123620148674334;
% x0=[1.000000000000000;0;0.056667402379709;0;0.591715153162376;0]; 2*1.310229401292405;
% 1.000 to 1.12 and 1.120 to 
% x0=[1.010000000000000;0;0.065589811196310;0;0.523743735947168;0]; 2*1.403831515057635;
% x0=[1.020000000000000;0;0.070976852313225;0;0.476362046334244;0]; 2*1.467552517073913;
% x0=[1.030000000000000;0;0.074102659855593;0;0.437901519523527;0]; 2*1.515414079199118;
% x0=[1.040000000000000;0;0.075473132397941;0;0.404214695876538;0]; 2*1.553278768721214;
% x0=[1.050000000000000;0;0.075315589841332;0;0.373341407980357;0]; 2*1.584246994771595;
% x0=[1.060000000000000;0;0.073712534600428;0;0.344187296764986;0]; 2*1.610181911905379;
% x0=[1.070000000000000;0;0.070641279606239;0;0.316064435304881;0]; 2*1.632294262170212;
% x0=[1.080000000000000;0;0.065968748259724;0;0.288495733099333;0]; 2*1.651411073694833;
% x0=[1.090000000000000;0;0.059396444769747;0;0.261118129557124;0]; 2*1.668115021372033;
% x0=[1.100000000000000;0;0.050283987170671;0;0.233627348644071;0]; 2*1.682822973060651;
% x0=[1.110000000000000;0;0.033430367143147;0;0.208995599183687;0];
% x0=[1.120000000000000;0;0.007307826276474;0;0.177162452569177;0]; 2*1.707348612623153;
% x0=[1.130000000000000;0;-0.176671768428050;0;-0.225470705431805;0]; 2*1.505051050732408;
% x0=[1.140000000000000;0;-0.163379713022233;0;-0.223431529531332;0]; 2*1.554095305416990;
% x0=[1.150000000000000;0;-0.146475086221399;0;-0.218155504051273;0]; 2*1.596818509505911;
% x0=[1.160000000000000;0;-0.124716788503581;0;-0.208673459534118;0]; 2*1.634496137789988;
% x0=[1.170000000000000;0;-0.103282995552828;0;-0.219409240867745;0]; 3.282456130305564;
% x0=[1.180000000000000;0;-0.029864233822740;0;-0.160820281841361;0]; 2*1.704107931530056;
% x0=[-1.600000000000000;0;0.638606243459622;0;1.205554648731255;0]; 2*3.118727703668909;
% x0=[-1.590000000000000;0;0.670017433586689;0;1.197923546134239;0]; 2*3.118631317018376;
% x0=[-1.580000000000000;0;0.699878909885960;0;1.190295963975811;0]; 2*3.118533135044251
% x0=[-1.005062644785412;0;1.574236144474778;0;0.756559543838613;0]; 2*3.108380902042432;
% -1.58 to -1.51 to -0.17
% x0=[-0.100000000000000;0;1.878718010024969;0;0.072172364279389;0]; 2*2.876829925605418;
%------------------------------- Needs to check ---------------------------
% {x0 = [0.834003975390563;-0.022000000000000;0;0;0;
%   0.000000000000000];
%------------------------------- Needs to check ---------------------------

% Used differential correction
[x0po, T] = poFam3BP3d(mu, eqNum, Ax1, Ax2, nFam);

%RelTol = 3.e-14 ; AbsTol = 1.e-16; % high accuracy;
%RelTol = 3.e-10 ; AbsTol = 1.e-12; % low accuracy
RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy

OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

% if saved data exist
% load x0po_T_L1.dat
% x0po = x0po_T_L1(:,1:6);
% T = x0po_T_L1(:,7);

figure()
hold on
for k = nFam:-1:1
    x0 = x0po(k,:) ;
%     x0=[-0.030000000000000;0;1.631848651200294;0;0.019706176807487;0];
    tf = T(k) ;
%     tf = 2*4.646062164537867/2;
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    plot3(x(:,1),x(:,2),x(:,3),'m.-', 'MarkerSize', 2);
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    hold on
    pause(0.01) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on
% plot3(0.836915146106164,0,0,'bs')
% plot3(1.155733835140644,0,0,'bs')
% plot3(-1.005062644785412,0,0,'bs')
% plot3(1-mu,0,0,'ro')
% plot3(-mu,0,0,'ro')
title("Halo orbits Differential Correction")

%% uncomment some of them to run the code
load x0po_T_L1_1.dat
load x0po_T_L1_10.dat
load x0po_T_L1_2.dat
load x0po_T_L1_3.dat
load x0po_T_L1_4.dat
load x0po_T_L1_5.dat
load x0po_T_L1_6.dat
load x0po_T_L1_7.dat
load x0po_T_L1_8.dat
load x0po_T_L1_9.dat
x0po = [x0po_T_L1_1; x0po_T_L1_2; x0po_T_L1_3; x0po_T_L1_4;...
    x0po_T_L1_5; x0po_T_L1_6; x0po_T_L1_7; x0po_T_L1_8;...
    x0po_T_L1_9; x0po_T_L1_10]; % 1-50, 51-end

load x0po_T_L2_1.dat
load x0po_T_L2_2.dat
load x0po_T_L2_3.dat
load x0po_T_L2_4.dat
load x0po_T_L2_5.dat
load x0po_T_L2_6.dat
load x0po_T_L2_7.dat
load x0po_T_L2_8.dat
load x0po_T_L2_9.dat
x0po = [x0po_T_L2_1; x0po_T_L2_2; x0po_T_L2_3; x0po_T_L2_4;...
    x0po_T_L2_5; x0po_T_L2_6; x0po_T_L2_7; x0po_T_L2_8;...
    x0po_T_L2_9];

load x0po_T_L3_13.dat
load x0po_T_L3.dat
load x0po_T_L3_2.dat
load x0po_T_L3_3.dat
load x0po_T_L3_4.dat
load x0po_T_L3_5.dat
load x0po_T_L3_6.dat
load x0po_T_L3_7.dat
load x0po_T_L3_8.dat
load x0po_T_L3_9.dat
load x0po_T_L3_10.dat
load x0po_T_L3_11.dat
load x0po_T_L3_12.dat
x0po = [x0po_T_L3; x0po_T_L3_2; x0po_T_L3_3;...
    x0po_T_L3_4; x0po_T_L3_5; x0po_T_L3_6;...
    x0po_T_L3_7; x0po_T_L3_8; x0po_T_L3_9;...
    x0po_T_L3_10; x0po_T_L3_11; x0po_T_L3_12;...
    x0po_T_L3_13];


figure()
hold on
% plot3(0.836915146106164,0,0,'ks')
% plot3(1.155733835140644,0,0,'ks')
plot3(-1.005062644785412,0,0,'ks')
% plot3(1-mu,0,0,'ko')
plot3(-mu,0,0,'ko')
grid on
for k = 1:10:size(x0po,1)
    x0 = x0po(k,1:6) ;
    tf = x0po(k, end) ;
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    plot3(x(:,1),x(:,2),x(:,3),'b.-', 'MarkerSize', 2);
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    hold on
    pause(0.01) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
title("Halo orbits, Differential Correction, L3")

%% initial condition: L1 with different velocity
figure()
hold on
grid on
for k = 1
    x0 = [0.836915146106164,0,0,0,0,k];
    tf = 400 ;
    [x,t] = trajGet3BP3d(x0,tf,mu,OPTIONS) ;
    plot3(x(:,1),x(:,2),x(:,3),'b.-', 'MarkerSize', 2);
    plot3(x(1,1),x(1,2),x(1,3),'r*', 'MarkerSize', 2);
    plot3(x(end,1),x(end,2),x(end,3),'ro', 'MarkerSize', 2);
    pause(0.01) ;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3(1-mu,0,0,'ko')
plot3(-mu,0,0,'ko')