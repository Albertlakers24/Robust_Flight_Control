%% Part #1 System Modeling
clc;
clearvars;

g = 9.80665; % m/s^2 Standard gravity acceleration
a = 316.0561; % m/s Speed of sound at 6096m
M = 3; % Mach number
Z_alpha = 1236.8918; % m/s^2 normal force derivative 
M_alpha = -300.4211; % 1/s^2 pitch moment derivative
M_q = 0; % 1/s demping derivative
Z_delta = 108.1144; % m/s^2 Control force derivative
M_delta = -131.3944; % 1/s^2 control moment derivative
A_alpha = 1434.7783; % m/s^2/rad normal acceleration derivative
A_delta = 115.0529; % m/s^2/rad control acceleration derivative 
omega_a = 150; % rad/s actuator natural frequency
zeta_a = 0.7; % actuator damping
V = M * a; % Airspeed
s = tf('s');

% Missile model
A_m = [-Z_alpha/V 1; M_alpha M_q]; % A matrix missile
B_m = [-Z_delta/V; M_delta;]; % B Matrix missile
C_m = [-A_alpha/g 0; 0 1]; % C Matrix missile
D_m = [-A_delta/g; 0]; % D matrix missile

G_m = ss(A_m, B_m, C_m, D_m);
G_m.InputName = 'u_m';
G_m.StateName = {'x1', 'x2'};
G_m.OutputName = {'y1', 'y2'};

% Actuator model
A_a = [0 1; -omega_a^2 -2*zeta_a*omega_a]; % A matrix actuator
B_a = [0; omega_a^2;]; % B Matrix actuator
C_a = [1 0; 0 1]; % C Matrix actuator
D_a = [0; 0]; % D matrix actuator

G_a = ss(A_a, B_a, C_a, D_a);
G_a.InputName = 'u_cmd';
G_a.StateName = {'x3', 'x4'};
G_a.OutputName = {'u_m', 'udot_m'};

% Save the models
save G_a G_a;
save G_m G_m;
load('G_a.mat');
load('G_m.mat');

linsys = linearize('Airframe');
T = [0 0 1 0;
     0 0 0 1;
     1 0 0 0;
     0 1 0 0];
G_am = ss2ss(linsys,T);

G_am_y1_u_cmd = G_am('y1', 'u_cmd');
G_am_y2_u_cmd = G_am('y2', 'u_cmd');

zpk_y1_u_cmd = zpk(G_am_y1_u_cmd);
zpk_y2_u_cmd = zpk(G_am_y2_u_cmd);

% figure;
% iopzmap(G_am);
% grid on;
% title('Pole-Zero Map');

%% Part #2a Loop Shaping
% % Damping gain tuning
% figure;
% rlocusplot(-G_am(2,1));
% sgrid(0.7,[]);
% title('Root Locus Plot');

C_q_gain = tf(-0.163); % C_q gain obtained from root locus plot for CL damping of 0.7
linsys_cq = linearize('ClosedLoop_Cq');
G_cl_q_unsc = ss2ss(linsys_cq,T);
zpk_G_cl_q_unsc = zpk(G_cl_q_unsc);

%% Part #2b Loop Shaping
% Scaling gain Design
C_sc_gain = tf(1/ dcgain(G_cl_q_unsc));
linsys_cq_csc = linearize("ClosedLoop_CqCsc");
G = ss2ss(linsys_cq_csc,T);
zpk_G = zpk(G);

% figure;
% step(G);
% grid on;
% title("Step Reponse of G");

%% Part #2c Loop Shaping
% Integral gain Design
C_i_gain = tf(1); % momentary unitary value
linsys_cq_csc_ci = linearize("ClosedLoop_CqCscCi");
T_cq_csc_ci = [0 0 0 1 0;
     0 0 0 0 1;
     0 1 0 0 0;
     0 0 1 0 0;
     1 0 0 0 0];
G_ol_nz = ss2ss(linsys_cq_csc_ci,T_cq_csc_ci);
%sisotool(G_ol_nz);
C_i_gain_phase = 6.2693;
C_i_gain_settling = 5.9043;
C_i_gain = tf(min(C_i_gain_settling,C_i_gain_phase));

%% Part #3a Weighting filters
%Weight 1 information
dcgain_1 = db2mag(-60);
freq_1 = 4;
mag_1 = db2mag(-3.01);
PM = pi/6;
M_1 = 1/(2*sin(PM/2));
hfgain_1 = M_1;

%Weight 2 information
dcgain_2 = db2mag(60);
freq_2 = omega_a;
mag_2 = db2mag(-15);
hfgain_2 = dcgain_1;

W_1_inv = makeweight(dcgain_1, [freq_1,mag_1], hfgain_1); 
W_2_inv = makeweight(dcgain_2, [freq_2,mag_2], hfgain_2);

W_1 = inv(W_1_inv);
W_2 = inv(W_2_inv);

%Recomputation of A1, A2, M1, M2, omega1 and omega2
[numW1, denW1] = tfdata(W_1,'v');
M1 = 1/numW1(1);
w1 = numW1(2);
A1 = denW1(2)/w1;

[numW2, denW2] = tfdata(W_2,'v');
M2 = denW2(1);
w2 = denW2(2);
A2 = w2/numW2(2);

% figure;
% bode(W_1, W_2)
% legend('W1','W2');
% grid on;

%% Part #3b Reference model
% Required constants
z_m         = zero(G);
z_m         = z_m(2);           % Non-minimum phase zero

%settlingtime_des= 0.18;             % 5% settle time [sec]
%overshoot_des  = 0.05;             % Overshoot

% Defining ranges for omega_d and zeta_d and vars to store results
%omega_d_range     = linspace(0.01, 100, 100);
%zeta_d_range = linspace(0.1, 1, 100);

%omega_d_refm  = 0;
%zeta_d_refm = 0;
%error       = inf;

%Brute force to find required omega_d and zeta_d
%for omega_d = omega_d_range
%    for zeta_d = zeta_d_range
%        num = omega_d^2 * [-1/z_m, 1];
%        den = [1, 2*zeta_d*omega_d, omega_d^2];
%        T_d = tf(num, den);
%        
%        info                = stepinfo(T_d, 'SettlingTimeThreshold', 0.05);
%        settling_time_error = abs(info.SettlingTime - settlingtime_des);
%        overshoot_error     = abs(info.Overshoot/100 - overshoot_des);
%
%        total_error         = settling_time_error + overshoot_error;
%        
%        if total_error < error
%            error = total_error;
%            omega_d_refm = omega_d;
%            zeta_d_refm = zeta_d;
%        end
%    end
%end


% Updating reference model with best parameters and plotting step response
%num = omega_d_refm^2 * [-1/z_m, 1];
%den = [1, 2*zeta_d_refm*omega_d_refm, omega_d_refm^2];
%T_d = zpk(tf(num, den)); 

%From the brute force method we get the next zeta_d, omega_d and TF
omega_d_refm = 18.1900;
zeta_d_refm = 0.7000;
num = omega_d_refm^2 * [-1/z_m, 1];
den = [1, 2*zeta_d_refm*omega_d_refm, omega_d_refm^2];
T_d = zpk(tf(num,den));

%% Part #3c.1 Controller design (hinfsyn case)
clc;
W_3 = W_1;
P = linearize("Design");
T_P = [0 0 1 0 0 0 0 0 0;
       0 0 0 1 0 0 0 0 0;
       1 0 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0 0;
       0 0 0 0 1 0 0 0 0;
       0 0 0 0 0 1 0 0 0;
       0 0 0 0 0 0 1 0 0;
       0 0 0 0 0 0 0 1 0;
       0 0 0 0 0 0 0 0 1];
P = ss2ss(P,T_P);
n_meas = 1;
n_cont = 1;
hinf_options = hinfsynOptions("RelTol",1e-6);
[C_e,T_wz,gamma] = hinfsyn(P,n_meas,n_cont,hinf_options);
sigma_options = sigmaoptions('cstprefs');
sigma_options.MagUnit = 'abs';
sigma_options.XLim = [1e-3, 1e4];
T_wz_inf = norm(T_wz,"inf");
if T_wz_inf <= gamma
    disp("The Maximum Singular Value is smaller than Global Performance Level")
end

% figure;
% sigma(zpk(T_wz),sigma_options);
% yline(gamma,'r--');
% yline(T_wz_inf,'b--');
% hold on;
% sigma(zpk(T_wz(1)),sigma_options);
% hold on;
% sigma(zpk(T_wz(2)),sigma_options);
% hold on;
% sigma(zpk(T_wz(3)),sigma_options);
% legend('Global T_{wz}', "Global \gamma","T_{wz} inf","T_{wz1}", "T_{wz2}", "T_{wz3}")

% Weight 3 information Re-evluation
dcgain_3 = db2mag(-60);
freq_3 = 4;
hfgain_3 = M_1;
% i = 0;
% while gamma < 1
%     mag_3 = db2mag(-22.641-i)
%     i = i + 0.0001
%     W_3_inv = makeweight(dcgain_3, [freq_3,mag_3], hfgain_3); 
%     W_3 = inv(W_3_inv);
% 
%     % P matrix reconstruction
%     P = ss2ss(linearize("Design"),T_P);
%     [C_e_reeval,T_wz_reeval,gamma_reeval] = hinfsyn(P,n_meas,n_cont,hinf_options);
%     gamma_1 = norm(T_wz_reeval(1),'inf');
%     gamma_2 = norm(T_wz_reeval(2),'inf');
%     gamma_3 = norm(T_wz_reeval(3),'inf');
%     gamma = max([gamma_1,gamma_2,gamma_3]);
% end
mag_3 = db2mag(-22.6415);
W_3_inv = makeweight(dcgain_3, [freq_3,mag_3], hfgain_3);
W_3 = inv(W_3_inv);
P = linearize("Design");
[C_e_reeval,T_wz_reeval,gamma_reeval] = hinfsyn(P,n_meas,n_cont,hinf_options);
gamma_1 = norm(T_wz_reeval(1),'inf');
gamma_2 = norm(T_wz_reeval(2),'inf');
gamma_3 = norm(T_wz_reeval(3),'inf');
T_wz_inf_reeval = norm(T_wz_reeval,"inf");
% figure;
% sigma(zpk(T_wz_reeval),sigma_options);
% yline(gamma_reeval,'r--');
% yline(T_wz_inf_reeval,'b--');
% hold on;
% sigma(zpk(T_wz_reeval(1)),sigma_options);
% hold on;
% sigma(zpk(T_wz_reeval(2)),sigma_options);
% hold on;
% sigma(zpk(T_wz_reeval(3)),sigma_options);
% grid on;
% legend('Global T_{wzr}', "Global \gamma","T_{wz} inf","T_{wz1}", "T_{wz2}", "T_{wz3}")
% title("Re-Evaluated Singular Value")

%% #3c.2 Controller order reduction
C0_e = C_e_reeval;
[z,p,k] = zpkdata(C0_e,'v');
k_min = k * z(1)/ p(1);
C_e_min = zpk(z(2:7),p(2:8),k_min);

% figure;
% bode(C0_e)
% hold on;
% bode(C_e_min);
% legend("Original System", "Minimized System");

C_i_min = zpk(C_e_min.Z{1}, C_e_min.P{1}(1:6), k_min);
C_i_red = balred(C_i_min,2);

% figure;
% bode(C_i_min);
% hold on
% bode(C_i_red);
% legend("C_i_min", "C_i_red");

%% Part #3c.3 Controller analysis & simulation
F_f = tf(1);
C_e_red = C_i_red * 1/s;
C_i_gain = C_i_red;
T = ss2ss(linearize("ClosedLoop_Test"),T_P);
So = zpk(T(1,1));
Ce_So = zpk(T(2,1));
To = zpk(T(3,1));
Tm = zpk(T(4,1));
Ti = zpk(-T(2,2));
SoG = zpk(T(3,2));
Si = zpk(T(5,2));

% figure;
% subplot(2,3,1);
% sigma(So,sigma_options,'b');
% hold on;
% sigma(Si,sigma_options,'b');
% hold on;
% sigma(W_1_inv,'r');
% grid on;
% title("S_{O}, S_{i}, W_{1}^{-1} Singular Value Plot");
% 
% subplot(2,3,2);
% sigma(Ce_So,sigma_options,'b');
% hold on;
% sigma(W_2_inv,'r');
% grid on;
% title("C_{e}S_{O}, W_{2}^{-1} Singular Value Plot");
% 
% subplot(2,3,3);
% sigma(Tm,sigma_options,'b');
% hold on;
% sigma(W_3_inv,'r');
% grid on;
% title("T_{m}, W_{3}^{-1} Singular Value Plot");
% 
% subplot(2,3,4);
% sigma(Ti,sigma_options,'b');
% hold on;
% sigma(To,'b');
% grid on;
% title("T_{i}, T_{O} Singular Value Plot");
% 
% subplot(2,3,5);
% sigma(SoG,sigma_options,'b');
% grid on;
% title("S_{O}G Singular Value Plot");
% 
% subplot(2,3,6);
% sigma(C_e,sigma_options,'b');
% hold on;
% sigma(C_e_red,sigma_options,'b');
% grid on;
% title("C_{e}, C_{e}_{red} Singular Value Plot");

T_OLT = [0 0 0 0 0 1 0;
         0 0 0 0 0 0 1;
         0 0 0 1 0 0 0;
         0 0 0 0 1 0 0;
         1 0 0 0 0 0 0;
         0 1 0 0 0 0 0;
         0 0 1 0 0 0 0];
sys_OLT = ss2ss(linearize("OpenLoop_Test"),T_OLT);
[Gm,Pm,Wcg,Wcp] = margin(sys_OLT);
Dm = deg2rad(Pm)/Wcp;

% figure;
% bode(sys_OLT);
% 
% figure;
% grid on;
% subplot(2,2,1);
% step(So,'b');
% 
% subplot(2,2,2);
% step(T_d,'r');
% hold on;
% step(To,'b');
% 
% subplot(2,2,3);
% step(SoG,'b');
% 
% subplot(2,2,4);
% step(zpk(T(6,1)),'b');

%% Part #3d Feedback controller desgin (hinfstruct case)
%3D.1
s = tf('s');
P = linearize('Design');
C_e_red_D0 = tunableTF('C_e_red_D0',2,2)*1/s;
hinfstruct_options = hinfstructOptions('UseParallel',true,'TolGain',10e-6);
[C_e_red_star, gamma_star, info_star] = hinfstruct(P, C_e_red_D0, hinfstruct_options);

C_e_red_star = zpk(C_e_red_star);
C_i_red_star = zpk(C_e_red_star)*s;

%Q1 5x2 table 
data = [27.01, 26.905; 
        36.32, 28.8; 
        689.7, 466.7; 
        62.38, 60.07; 
        2258, 1526];
rowNames = {'K', 'n_1', 'n_2', 'd_1', 'd_2'};
columnNames = {'C_i_red', 'C_i_red_star'};
Table5x2 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

%Q2
T_wz_D = lft(P,C_e,1,1);
T_wz_D_1 = tf(T_wz_D(1));
T_wz_D_2 = tf(T_wz_D(2));
T_wz_D_3 = tf(T_wz_D(3));

gamma_D_1 = norm(T_wz_D_1,'inf');
gamma_D_2 = norm(T_wz_D_2,'inf');
gamma_D_3 = norm(T_wz_D_3,'inf');

figure(), hold on
sigma(T_wz_D,sigma_options)
sigma(T_wz_D_1,sigma_options)
sigma(T_wz_D_2,sigma_options)
sigma(T_wz_D_3,sigma_options)
grid on
legend('T_{wz}','T_{wz(1)}','T_{wz(2)}','T_{wz(3)}')

%plot but in DB
%figure(), hold on
%sigma(T_wz_D)
%sigma(T_wz_D_1)
%sigma(T_wz_D_2)
%sigma(T_wz_D_3)
%grid on
%legend('T_{wz}','T_{wz(1)}','T_{wz(2)}','T_{wz(3)}')

figure(), hold on
bode(C_i_min,'b')
bode(C_i_red,'g')
bode(C_i_red_star,'m')
legend('C-i-min','C-i-red','C-i-red-star')

%3D.2
C_i_gain = C_i_red_star;
T_3D = linearize("ClosedLoop_Test");

So_3D = zpk(T_3D(1,1));
Ce_So_3D = zpk(T_3D(2,1));
To_3D = zpk(T_3D(3,1));
Tm_3D = zpk(T_3D(4,1));
Ti_3D = zpk(-T_3D(2,2));
SoG_3D = zpk(T_3D(3,2));
Si_3D = zpk(T_3D(5,2));


figure;
subplot(2,3,1);
sigma(So_3D,sigma_options,'m');
hold on;
sigma(Si_3D,sigma_options,'m');
hold on;
sigma(So,sigma_options,'b');
hold on;
sigma(Si,sigma_options,'b');
hold on;
sigma(W_1_inv,'r');
grid on;
title("Singular Value Plot 3D");
legend('S_{O,hinfstruct}','S_{i,hinfstruct}','S_{O,hinfsyn}','S_{i,hinfsyn}','W_{1}^{-1}')
hold off;

subplot(2,3,2);
sigma(Ce_So_3D,sigma_options,'m');
hold on;
sigma(Ce_So,sigma_options,'b');
hold on;
sigma(W_2_inv,'r');
grid on;
title("Singular Value Plot 3D");
legend('C_{e}S_{0,hinfstruct}','C_{e}S_{0,hinfsyn}','W_{2}^{-1}')
hold off;

subplot(2,3,3);
sigma(Tm_3D,sigma_options,'m');
hold on;
sigma(Tm,sigma_options,'b');
hold on;
sigma(W_3_inv,'r');
grid on;
title("Singular Value Plot 3D");
legend('T_{m,hinfstruct}','T_{m,hinfsyn}','W_{3}^{-1}')
hold off;

subplot(2,3,4);
sigma(Ti_3D,sigma_options,'m');
hold on;
sigma(To_3D,sigma_options,'m');
hold on;
sigma(Ti,sigma_options,'b');
hold on;
sigma(To,sigma_options,'b');
grid on;
title("Singular Value Plot 3D");
legend('T_{i,hinfstruct}','T_{O,hinfstruct}','T_{i,hinfsyn}','T_{o,hinfsyn}','W_{1}^{-1}')
hold off;


subplot(2,3,5);
sigma(SoG_3D,sigma_options,'m');
hold on;
sigma(SoG,sigma_options,'b');
grid on;
title("Singular Value Plot 3D");
legend('S_o_{G,hinfstruct}','S_o_{G,hinfsyn}')
hold off;


subplot(2,3,6);
sigma(C_e_red_star,sigma_options,'m');
hold on;
sigma(C_e_red,sigma_options,'b');
hold on;
sigma(C_e,sigma_options,'r');
grid on;
title("Singular Value Plot 3D");
legend('C_e_red_star','C_e_red','C_e')
hold off;

%second exercise
sys_OLT_D = linearize("OpenLoop_Test");
[Gm_D,Pm_D,Wcg_D,Wcp_D] = margin(sys_OLT_D);
Dm = deg2rad(Pm_D)/Wcp_D;

figure;
bode(sys_OLT_D,'m');
hold on;
bode(sys_OLT,'b');

figure;
grid on;
subplot(2,2,1);
step(So_3D,'m');
hold on;
step(So,'b');
 
subplot(2,2,2);
step(T_d,'r');
hold on;
step(To_3D,'m');
hold on;
step(To,'b');

 
subplot(2,2,3);
step(SoG_3D,'m');
hold on;
step(SoG,'b');
 
subplot(2,2,4);
step(zpk(T_3D(6,1)),'m');
hold on;
step(zpk(T(6,1)),'b');

%% Part #3eFeedforward controller design
F_f_init = zpk(T_d * minreal(To_3D)^(-1));

sigma_options.XLim = [1e-3, 1e4];
sigma_options.MagUnits = 'abs';

figure;
sigma(F_f_init,'b');

[z,p,k] = zpkdata(F_f_init,'v');
z([1,3,4]) = [];
p([3]) = [];
F_f_lf = zpk(z,p,k);
dcgain_init = dcgain(F_f_init);
dcgain_lf = dcgain(F_f_lf);
adjusted_k = k * dcgain_init / dcgain_lf;
F_f_lf = zpk(z,p,adjusted_k);

Freq_interval = [0.05, 100];
opts = balredOptions('StateElimMethod' , 'Truncate', "FreqIntervals", Freq_interval, "ErrorBound", "abs");
F_f = zpk(balred(F_f_lf,2,opts));

figure;
subplot(3,1,1)
sigma(F_f_init,'b',sigma_options);
title("Singular Value of F_f_init")

subplot(3,1,2)
sigma(F_f_lf,"r",sigma_options);
title("Singular Value of F_f_lf")

subplot(3,1,3)
sigma(F_f,"g",sigma_options);
title("Singular Value of F_f")



%% Part #4 Feedback controller redesign (systune case) (Optional)

