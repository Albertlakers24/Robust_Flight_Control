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


figure;
iopzmap(G_am);
grid on;
title('Pole-Zero Map');

%% Part #2a Loop Shaping
% Damping gain tuning
clc;
% figure;
% rlocusplot(-G_am(2,1));
% title('Root Locus Plot');
% sgrid(0.7,[]);
C_q_gain = -0.163; % C_q gain obtained from root locus plot for CL damping of 0.7
linsys_cq = linearize('ClosedLoop_Cq');
G_cl_q = ss2ss(linsys_cq,T);
G_cl_q_unsc = G_cl_q('y1','u_unsc');
zpk_G_cl_q_unsc = zpk(G_cl_q_unsc);

%% Part #2b Loop Shaping
% Scaling gain Design
C_sc_gain = 1/ dcgain(G_cl_q_unsc);
linsys_cq_csc = linearize("ClosedLoop_CqCsc");
G = ss2ss(linsys_cq_csc,T);
G = G('y1','u_p');
zpk_G = zpk(G);
G_tf = tf(G);
% figure;
% step(G_tf);
% grid on;
% title("Step Reponse of G");

%% Part #2c Loop Shaping
% Integral gain Design
C_i_gain = 1; % momentary unitary value
linsys_cq_csc_ci = linearize("ClosedLoop_CqCscCi");
T_cq_csc_ci = [0 0 1 0 0;
     0 0 0 1 0;
     1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 0 1];
G_ol_nz = ss2ss(linsys_cq_csc_ci,T_cq_csc_ci);
% sisotool(G_ol_nz);
C_i_gain_phase = 6.2693;
C_i_gain_settling = 5.4738;
C_i_gain = min(C_i_gain_settling,C_i_gain_phase);

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
T_d = zpk(tf(num, den));


%% Part #3c Feedback controller design (hinfsyn case)
clc;
W_3 = W_1;
P = [W_1 -W_1*zpk_G;
    0 W_2;
    W_3*T_d -W_3*zpk_G;
    1 -zpk_G];
n_meas = 1;
n_cont = 1;
hinf_options = hinfsynOptions;
hinf_options.RelTol = 1e-6;
[C_e,T_wz,gamma] = hinfsyn(P,n_meas,n_cont,hinf_options);
sigma_options = sigmaoptions('cstprefs');
sigma_options.MagUnit = 'abs';
sigma_options.YLim = [-0.2,1.2];
[sv, freq] = sigma(zpk(T_wz));
T_wz_inf = norm(T_wz,"inf");
if T_wz_inf <= gamma
    disp("The Maximum Singular Value is smaller than Global Performance Level")
end

figure;
sigma(zpk(T_wz),sigma_options);
yline(gamma,'r--');
yline(T_wz_inf,'b--');
hold on;
sigma(zpk(T_wz(1)),sigma_options);
hold on;
sigma(zpk(T_wz(2)),sigma_options);
hold on;
sigma(zpk(T_wz(3)),sigma_options);
legend('Global T_{wz}', "Global \gamma","T_{wz} inf","T_{wz1}", "T_{wz2}", "T_{wz3}")

%Weight 3 information Re-evluation
dcgain_3 = db2mag(-60);
freq_3 = 4;
hfgain_3 = M_1;
i = 0;
while gamma < 1
    mag_3 = db2mag(-22.64-i);
    W_3_inv = makeweight(dcgain_3, [freq_3,mag_3], hfgain_3); 
    W_3 = inv(W_3_inv);

    %P matrix reconstruction
    P = [W_1 -W_1*zpk_G;
        0 W_2;
        W_3*T_d -W_3*zpk_G;
        1 -zpk_G];
    [C_e_reeval,T_wz_reeval,gamma_reeval,info_reeval] = hinfsyn(P,n_meas,n_cont,tolerance);
    gamma_1 = norm(T_wz_reeval(1),'inf');
    gamma_2 = norm(T_wz_reeval(2),'inf');
    gamma_3 = norm(T_wz_reeval(3),'inf');
    gamma = max([gamma_1,gamma_2,gamma_3]);
    i = i + 0.0001;
end
% mag_3 = -22.6409
figure;
sigma(zpk(T_wz_reeval),sigma_options);
yline(gamma_reeval,'r--');
yline(T_wz_inf_reeval,'b--');
hold on;
sigma(zpk(T_wz_reeval(1)),sigma_options);
hold on;
sigma(zpk(T_wz_reeval(2)),sigma_options);
hold on;
sigma(zpk(T_wz_reeval(3)),sigma_options);
legend('Global T_{wzr}', "Global \gamma","T_{wz} inf","T_{wz1}", "T_{wz2}", "T_{wz3}")
title("Re-evaluated Singular Value")

%% Part #3d Feedback controller desgin (hinfstruct case)

%% Part #3eFeedforward controller design

%% Part #4 Feedback controller redesign (systune case) (Optional)
