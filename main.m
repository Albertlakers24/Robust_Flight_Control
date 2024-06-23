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
linsys_csc = linearize("ClosedLoop_CqCsc");
G = ss2ss(linsys_csc,T);
G = G('y1','u_p');
zpk_G = zpk(G);
G_tf = tf(G);
figure;
step(G_tf);
grid on;
title("Step Reponse of G");

%% Part #3a Weighting filters

%% Part #3b Reference model

%% Part #3c Feedback controller design (hinfsyn case)

%% Part #3d Feedback controller desgin (hinfstruct case)

%% Part #3eFeedforward controller design

%% Part #4 Feedback controller redesign (systune case) (Optional)