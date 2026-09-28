%% Parameters (ideal SEPIC)
clear;
clc;
P = 10;
Vo = 20;
Vin = 12;
L1 = 300e-6;
L2 = 500e-6;
Cs = 5.2e-6;
Co = 15.625e-6;
Io = P/Vo;              % 0.5 A
R = Vo/Io;               % 40 ohm
d = (Vo/Vin)/(1 + (Vo/Vin));   % 0.625
Vcs = Vin;                % 12 V 
il1 = P/Vin;              % 0.8333 A  
il2 = Io;                 % 0.5 A
Ts = 10e-6;
Vf = 0.3;
%% Constants used by the Simulink model
x_nom   = [il1; il2; Vcs; Vo];      % linearization point of the plant model
D_nom   = d;
Vin_nom = Vin;
x0_dev  = [0.1; -0.1; 0; 0];        % initial kick along the internal-mode shape
Vo_ref  = 20;   Po_ref = 10;
dmin    = 0.05; dmax = 0.95;
alpha   = 0.01;                     % Vin filter coefficient

%% small signal model
A_hat = [0        0       -(1-d)/L1  -(1-d)/L1;
         0        0        d/L2      -(1-d)/L2;
         (1-d)/Cs -d/Cs    0          0;
         (1-d)/Co  (1-d)/Co 0        -1/(R*Co)];

B_hat = [(Vcs+Vo)/L1   1/L1;
         (Vcs+Vo)/L2   0;
         -(il1+il2)/Cs 0;
         -(il1+il2)/Co 0];
%% Small-signal model with diode drop
A_hat = [0         0        -(1-d)/L1   -(1-d)/L1;
         0         0         d/L2       -(1-d)/L2;
         (1-d)/Cs  -d/Cs     0           0;
         (1-d)/Co  (1-d)/Co  0          -1/(R*Co)];

B_hat = [(Vcs+Vo+Vf)/L1   1/L1;
         (Vcs+Vo+Vf)/L2   0;
         -(il1+il2)/Cs    0;
         -(il1+il2)/Co    0];

C = [0 0 0 1];
D = [0 0];

sys = ss(A_hat, B_hat, C, D);
sysd_z = c2d(sys, Ts, 'zoh');
[Ad, Bd] = ssdata(sysd_z);
%% MPC

% select the control input for prediction

Bu = Bd(:,1);          % duty-cycle column only; Vin (column 2) is a disturbance, not a decision variable
nx = size(Ad,1);        % number of states = 4 for the ideal model
% Horizon
Np = 30;     % prediction horizon, in samples (30*10us = 300us)
Nu = Np;     % control horizon = Np for now; Nu < Np comes later

% Prediction matrices M and S
M = zeros(Np*nx, nx);
S = zeros(Np*nx, Nu);   % nu=1, so each column is one decision variable

Apow = eye(nx);
for i = 1:Np
    Apow = Ad * Apow;                      % Ad^i, built incrementally
    M((i-1)*nx+1 : i*nx, :) = Apow;        % stack Ad^1, Ad^2, ..., Ad^Np
end

for i = 1:Np
    for j = 1:i                            % only j <= i contributes (causality)
        S((i-1)*nx+1 : i*nx, j) = Ad^(i-j) * Bu;
    end
end

% Weights 
diL1_max = 0.10*il1;    % 83 mA  (10% of nominal)
diL2_max = 0.10*il2;    % 50 mA
dvCs_max = 0.5;         % 0.5 V
dvo_max  = 0.05;         % 0.05 V 

Q_bryson = diag(1 ./ [diL1_max, diL2_max, dvCs_max, dvo_max].^2);

% Step 4b - Internal-mode term: penalize (diL1 - diL2)^2 2(XtMtQSU) + 2(XtMtQ_modeSU) = 2(XtMt(Q + Q_mode)SU)
c_mode      = [1 -1 0 0];     % extracts diL1 - diL2
dmode_max   = 0.01;           % 50 mA tolerance on the imbalance
w_mode      = 1;              % 1 = on, 0 = off (for comparisons)
Q_mode      = w_mode * (c_mode' * c_mode) / dmode_max^2;

Q = Q_bryson + Q_mode;

% Step 4c - R: same idea, tolerance on duty-cycle deviation
dd_max = 0.025;                % 5% duty-cycle deviation
Rw     = 1 / dd_max^2;
% Step 5 - QP matrices (offline)
Qbar = kron(eye(Np), Q);
Rbar = kron(eye(Nu), Rw);

H = 2*(S'*Qbar*S + Rbar);
H = (H + H')/2;
Fgain = 2*S'*Qbar*M;          % online: f = Fgain*dx0