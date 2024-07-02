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

linsys = linearize('Airframenew');
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

%% Part #2 Loop Shaping

% Damping gain tuning
clc;
figure;
rlocusplot(-G_am(2,1));
title('Root Locus Plot');
sgrid(0.7,[]);
C_q_gain = -0.163; % C_q gain obtained from root locus plot for CL damping of 0.7
linsys_cq = linearize('ClosedLoop_Cq');
G_cl_q = ss2ss(linsys_cq,T);
G_cl_q_unsc = G_cl_q('y1','u_unsc');
zpk_G_cl_q_unsc = zpk(G_cl_q_unsc);

% Scaling gain Design
C_sc_gain = 1/ dcgain(G_cl_q_unsc);
linsys_cq_csc = linearize("ClosedLoop_CqCsc");
G = ss2ss(linsys_cq_csc,T);
G = G('y1','u_p');
zpk_G = zpk(G);
G_tf = tf(G);
figure;
step(G_tf);
grid on;
title("Step Reponse of G");

% Integral gain Design
C_i_gain = db2mag(20); % momentary unitary value
linsys_cq_csc_ci = linearize("ClosedLoop_CqCscCi");
T_cq_csc_ci = [0 0 1 0 0;
     0 0 0 1 0;
     1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 0 1];
G_ol_nz = ss2ss(linsys_cq_csc_ci,T_cq_csc_ci);
zpk_G_ol_nz = zpk(G_ol_nz);
figure;
bode(G_ol_nz);
grid on;
title("Bode Plot of G_ol_nz");
C_i_gain = 1; % momentary unitary value
linsys_cq_csc_ci = linearize("ClosedLoop_CqCscCi");
T_cq_csc_ci = [0 0 1 0 0;
     0 0 0 1 0;
     1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 0 1];
G_ol_nz = ss2ss(linsys_cq_csc_ci,T_cq_csc_ci);
sisotool(G_ol_nz);
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

bode(W_1, W_2)
legend('W1','W2');
grid on;





%% Part #3b Reference model
G_zeros = zero(G);
s = tf('s');
z_m = G_zeros(1);

M_d = 5;
ts = 0.18;
error = 0.1;
w_d_correct =NaN;
z_d_correct = NaN;

%for w_d = 18:0.01:22
%    for z_d = 0.65:0.0001:0.7
%        T_d = (w_d^2*(-s/z_m + 1))/(s^2 + 2*z_d*w_d*s + w_d^2);
%        if abs(stepinfo(T_d).Overshoot - max_os) <= error && abs(stepinfo(T_d).SettlingTime - ts) <= error
%             w_d_correct = w_d;
%             z_d_correct = z_d;
%        end
%    end
%end


%% Part #3c Feedback controller design (hinfsyn case)

%% Part #3d Feedback controller desgin (hinfstruct case)

%% Part #3eFeedforward controller design

%% Part #4 Feedback controller redesign (systune case) (Optional)
