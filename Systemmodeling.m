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
A_a = [0 1; -omega_a^2 -2*zeta_a*omega_a]; % A matrix missile
B_a = [0; omega_a^2;]; % B Matrix missile
C_a = [1 0; 0 1]; % C Matrix missile
D_a = [0; 0]; % D matrix missile

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