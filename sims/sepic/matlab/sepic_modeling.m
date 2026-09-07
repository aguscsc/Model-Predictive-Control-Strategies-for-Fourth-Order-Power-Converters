%% Parameters (ideal SEPIC)
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

%% small signal model
A_hat = [0        0       -(1-d)/L1  -(1-d)/L1;
         0        0        d/L2      -(1-d)/L2;
         (1-d)/Cs -d/Cs    0          0;
         (1-d)/Co  (1-d)/Co 0        -1/(R*Co)];

B_hat = [(Vcs+Vo)/L1   1/L1;
         (Vcs+Vo)/L2   0;
         -(il1+il2)/Cs 0;
         -(il1+il2)/Co 0];

C = [0 0 0 1];
D = [0 0];

sys = ss(A_hat, B_hat, C, D);

%% Tustin discretization
Ts = 10e-6;
sysd = c2d(sys, Ts, 'tustin');
[Ad, Bd, Cd, Dd] = ssdata(sysd);

%% optimization problem