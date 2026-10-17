%% MPC for the NON-IDEAL SEPIC
% Same MPC structure as sepic_mpc_sim.m (dense QP, deviation coordinates,
% Bryson weights + internal-mode term, adaptive reference), but the
% prediction model and the plant include the converter parasitics.
%
% States:  x = [iL1; iL2; vCs; vCo]   (vCs, vCo = voltages on the ideal
%          capacitance, i.e. behind the ESR; vo is an output, not a state)
% Inputs:  u = [Vin; Vf]               (Vf enters as a constant source)
% Current directions follow the ideal model: iL1 into the switch node,
% iL2 from ground up through L2 into the diode node.
% CCM is assumed throughout (both in the model and in the plant).
%
% Offset-free MPC (use_offset_free = true):
% The prediction model works in coordinates relative to the Vin-scheduled
% reference, e = x - x_ref(Vin), ue = d - D_ref(Vin), and is augmented
% with a constant disturbance w (one component per state):
%     e(k+1) = Ad*e(k) + Bu*ue(k) + w,      w(k+1) = w(k)
% Because all four states are measured, w is estimated from the one-step
% prediction error, low-pass filtered. Everything the model gets wrong
% (load changes, parasitic mismatch, switching ripple, linearisation)
% ends up in w_hat. A target calculation then finds the steady state
% (es, us) that gives vo = Vo despite w_hat, and the QP regulates to it.
clear;
clc;

%% ===================== PARAMETERS (edit here) ==========================
% --- Operating point ---
P    = 10;          % W, nominal output power
Vo   = 20;          % V, output voltage reference
Vin  = 12;          % V, nominal input voltage
R    = Vo^2/P;      % ohm, load (40 ohm)

% --- Reactive components ---
L1 = 300e-6;        % H
L2 = 500e-6;        % H
Cs = 5.2e-6;        % F, coupling capacitor
Co = 15.625e-6;     % F, output capacitor

% --- Parasitics ---
rL1    = 0.19;      % ohm, L1 winding DCR (~1.48 m of 0.404 mm Cu)
rL2    = 0.25;      % ohm, L2 winding DCR
Rsh1   = 0.05;      % ohm, current-sense shunt in series with L1 (0 if none)
Rsh2   = 0.05;      % ohm, current-sense shunt in series with L2 (0 if none)
rCs    = 0.01;      % ohm, ESR of Cs
rCo    = 0.01;      % ohm, ESR of Co
Rds_on = 17.5e-3;   % ohm, MOSFET on-resistance (IRFZ44N)
Vf     = 0.3;       % V,   diode forward-voltage (threshold part)
Rd     = 0.03;      % ohm, diode dynamic resistance

% --- Sampling ---
Ts = 10e-6;         % s, control period = switching period

% --- MPC horizon ---
Np = 30;            % prediction horizon (samples)
Nu = Np;            % control horizon (Nu <= Np; moves after Nu are held)

% --- MPC weights (Bryson tolerances) ---
diL1_max  = 0.10;   % fraction of nominal iL1
diL2_max  = 0.10;   % fraction of nominal iL2
dvCs_max  = 0.5;    % V
dvo_max   = 0.05;   % V
dmode_max = 0.01;   % A, tolerance on the internal-mode imbalance iL1 - iL2
w_mode    = 1;      % 1 = internal-mode term on, 0 = off
dd_max    = 0.025;  % duty-cycle deviation tolerance

% --- Offset-free (disturbance estimator + target calculation) ---
use_offset_free = true;   % false = plain MPC (offset under load changes)
beta_w          = 0.05;   % disturbance estimator filter coefficient per sample
                          % (time constant ~ Ts/beta_w; larger = faster, noisier)

% --- Constraints ---
dmin = 0.05;
dmax = 0.75;

% --- Reference generator ---
alpha      = 0.01;  % Vin low-pass filter coefficient
Vin_lut    = 8:0.05:16;  % Vin grid for the steady-state reference table

% --- Simulation ---
Tsim         = 10e-3;               % s
t_vin_step   = 3e-3;                % s, time of the Vin step
Vin_step     = 14;                % V, Vin after the step
t_load_step  = 6e-3;              % s, time of the load step (Inf = no step)
R_step       = 30;                  % ohm, load after the step
x0_kick      = [0; 0; 0; 0];   % initial deviation from the operating point
use_mpc      = true;                % false = open loop (d held at nominal
plant_model  = 'averaged';          % 'averaged' or 'switched' (on/off sub-intervals)
compare_ideal = false;               % Bode plot of ideal vs non-ideal d -> vo

%% Pack parameters
par = struct('L1',L1,'L2',L2,'Cs',Cs,'Co',Co,'R',R, ...
             'rL1',rL1+Rsh1,'rL2',rL2+Rsh2,'rCs',rCs,'rCo',rCo, ...
             'Rds',Rds_on,'Vf',Vf,'Rd',Rd);
[A1, B1, C1, A2, B2, C2] = sepic_switched_matrices(par);

%% Non-ideal operating point: find D that gives vo = Vo at Vin
[D_nom, x_nom] = sepic_operating_point(A1,B1,C1,A2,B2,C2, Vin, Vf, Vo);
il1 = x_nom(1);  il2 = x_nom(2);  vcs = x_nom(3);  vco = x_nom(4);

D_ideal = Vo/(Vo + Vin);
Pin  = Vin*il1;
Pout = Vo^2/R;
fprintf('Operating point (non-ideal):\n');
fprintf('  D = %.4f (ideal %.4f)\n', D_nom, D_ideal);
fprintf('  iL1 = %.4f A, iL2 = %.4f A, vCs = %.3f V, vCo = %.3f V\n', il1, il2, vcs, vco);
fprintf('  efficiency = %.1f %%\n', 100*Pout/Pin);

%% Small-signal model around (D_nom, x_nom)
U_nom = [Vin; Vf];
A_hat = D_nom*A1 + (1-D_nom)*A2;
Bd_c  = (A1 - A2)*x_nom + (B1 - B2)*U_nom;   % d -> x
Bv_c  = D_nom*B1(:,1) + (1-D_nom)*B2(:,1);   % Vin -> x
B_hat = [Bd_c Bv_c];

C_vo  = D_nom*C1 + (1-D_nom)*C2;             % vo = C(d)*x
D_vo  = (C1 - C2)*x_nom;                     % d -> vo feedthrough (ESR term)

sys_x  = ss(A_hat, B_hat, eye(4), zeros(4,2));
sysd_z = c2d(sys_x, Ts, 'zoh');
[Ad, Bd] = ssdata(sysd_z);

fprintf('Continuous poles (non-ideal):\n');  disp(eig(A_hat));

%% Ideal vs non-ideal control-to-output comparison
if compare_ideal
    d_i = D_ideal;  Ii1 = P/Vin;  Ii2 = P/Vo;
    A_id = [0         0        -(1-d_i)/L1  -(1-d_i)/L1;
            0         0         d_i/L2      -(1-d_i)/L2;
            (1-d_i)/Cs -d_i/Cs  0            0;
            (1-d_i)/Co (1-d_i)/Co 0         -1/(R*Co)];
    B_id = [(Vin+Vo)/L1; (Vin+Vo)/L2; -(Ii1+Ii2)/Cs; -(Ii1+Ii2)/Co];
    G_id = ss(A_id, B_id, [0 0 0 1], 0);
    G_ni = ss(A_hat, Bd_c, C_vo, D_vo);
    figure;
    bode(G_id, G_ni); grid on;
    legend('Ideal', 'Non-ideal');
    title('Control-to-output  d \rightarrow v_o');
end

%% MPC prediction matrices  X = M*dx(k) + S*U
Bu = Bd(:,1);      % duty cycle: decision variable
nx = size(Ad,1);

Apow = cell(Np+1, 1);
Apow{1} = eye(nx);
for n = 1:Np
    Apow{n+1} = Ad * Apow{n};
end

M = zeros(Np*nx, nx);
S = zeros(Np*nx, Nu);
for i = 1:Np
    rows = (i-1)*nx+1 : i*nx;
    M(rows, :) = Apow{i+1};
    for j = 1:i
        col = min(j, Nu);                         % moves beyond Nu are held at u(Nu)
        S(rows, col) = S(rows, col) + Apow{i-j+1} * Bu;
    end
end

%% Weights
Q_bryson = diag(1 ./ [diL1_max*il1, diL2_max*il2, dvCs_max, dvo_max].^2);
c_mode   = [1 -1 0 0];
Q_mode   = w_mode * (c_mode' * c_mode) / dmode_max^2;
Q        = Q_bryson + Q_mode;
Rw       = 1 / dd_max^2;

Qbar = kron(eye(Np), Q);
Rbar = kron(eye(Nu), Rw);

H = 2*(S'*Qbar*S + Rbar);
H = (H + H')/2;
Fgain = 2*S'*Qbar*M;          % online: f = Fgain*dx0

[~, p] = chol(H);
fprintf('H positive definite: %d | cond(H) = %.2e\n', p==0, cond(H));

%% Offset-free target calculation (offline)
% Solve for the steady state (es, us) of the model with disturbance w_hat
% that puts vo exactly on the reference:
%   [I - Ad, -Bu; C_vo, D_vo] * [es; us] = [w_hat; 0]
T_ss  = [eye(nx) - Ad, -Bu; C_vo, D_vo];
fprintf('Target-calculation matrix: cond = %.2e\n', cond(T_ss));
T_inv = inv(T_ss);
T_w   = T_inv(:, 1:nx);       % online: [es; us] = T_w * w_hat

%% Steady-state reference table vs Vin (non-ideal model, solved offline)
nL = numel(Vin_lut);
D_lut = zeros(1, nL);
X_lut = zeros(4, nL);
for n = 1:nL
    [D_lut(n), X_lut(:,n)] = sepic_operating_point(A1,B1,C1,A2,B2,C2, Vin_lut(n), Vf, Vo);
end

%% Closed-loop simulation
N = round(Tsim/Ts);
t = (0:N-1)*Ts;

vin_profile = Vin*ones(1,N);
vin_profile(t >= t_vin_step) = Vin_step;
R_profile = R*ones(1,N);
R_profile(t >= t_load_step) = R_step;

options = optimoptions('quadprog', 'Display', 'off');

x        = x_nom + x0_kick;      % physical (large-signal) plant state
vin_filt = Vin;
w_hat    = zeros(nx, 1);         % disturbance estimate
e_prev   = [];                   % previous e and ue, for the estimator
ue_prev  = [];

x_log    = zeros(4, N);
xref_log = zeros(4, N);
w_log    = zeros(4, N);
vo_log   = zeros(1, N);
d_log    = zeros(1, N);

for k = 1:N
    % --- measure ---
    vin_meas = vin_profile(k);

    % --- adaptive reference from the non-ideal steady-state table ---
    vin_filt = alpha*vin_meas + (1-alpha)*vin_filt;
    D_ref_k  = interp1(Vin_lut, D_lut, vin_filt, 'linear', 'extrap');
    x_ref_k  = interp1(Vin_lut, X_lut.', vin_filt, 'linear', 'extrap').';

    % --- disturbance estimator (one-step prediction error, filtered) ---
    e = x - x_ref_k;
    if use_offset_free && ~isempty(e_prev)
        w_raw = e - Ad*e_prev - Bu*ue_prev;
        w_hat = (1-beta_w)*w_hat + beta_w*w_raw;
    end

    % --- target calculation: shift the reference to cancel w_hat ---
    if use_offset_free
        tgt = T_w * w_hat;
        es  = tgt(1:nx);
        us  = tgt(nx+1);
    else
        es  = zeros(nx, 1);
        us  = 0;
    end
    x_tgt = x_ref_k + es;
    d_tgt = D_ref_k + us;

    % --- controller ---
    if use_mpc
        dx0 = x - x_tgt;
        f   = Fgain * dx0;
        lb  = (dmin - d_tgt) * ones(Nu,1);
        ub  = (dmax - d_tgt) * ones(Nu,1);
        [U_opt, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, [], options);
        if exitflag <= 0
            warning('QP failed at k = %d', k);
            U_opt = zeros(Nu,1);
        end
        d_k = d_tgt + U_opt(1);
    else
        d_k = D_nom;
    end
    e_prev  = e;
    ue_prev = d_k - D_ref_k;

    % --- plant (load may change; rebuild sub-circuit matrices if it did) ---
    if k == 1 || R_profile(k) ~= R_profile(k-1)
        par.R = R_profile(k);
        [A1p, B1p, C1p, A2p, B2p, C2p] = sepic_switched_matrices(par);
    end

    % --- log (vo uses the averaged output equation) ---
    x_log(:,k)    = x;
    xref_log(:,k) = x_tgt;
    w_log(:,k)    = w_hat;
    vo_log(k)     = (d_k*C1p + (1-d_k)*C2p) * x;
    d_log(k)      = d_k;

    % --- plant: advance one period with the large-signal non-ideal model ---
    u = [vin_meas; Vf];
    switch plant_model
        case 'averaged'
            x = lti_step(d_k*A1p + (1-d_k)*A2p, d_k*B1p + (1-d_k)*B2p, x, u, Ts);
        case 'switched'
            x = lti_step(A1p, B1p, x, u, d_k*Ts);       % switch on
            x = lti_step(A2p, B2p, x, u, (1-d_k)*Ts);   % switch off, diode on
        otherwise
            error('plant_model must be ''averaged'' or ''switched''');
    end
end

%% Plots
vo_ref_log = Vo*ones(1,N);
names = {'i_{L1} (A)', 'i_{L2} (A)', 'v_{Cs} (V)', 'v_o (V)'};
figure;
for n = 1:3
    subplot(5,1,n);
    plot(t*1e3, x_log(n,:), t*1e3, xref_log(n,:), '--'); grid on;
    ylabel(names{n});
end
subplot(5,1,4);
plot(t*1e3, vo_log, t*1e3, vo_ref_log, '--'); grid on;
ylabel(names{4});
subplot(5,1,5);
plot(t*1e3, d_log); grid on;
ylabel('d'); xlabel('t (ms)');
sgtitle(sprintf('Non-ideal SEPIC MPC (%s plant, offset-free = %d)', plant_model, use_offset_free));

if use_offset_free
    figure;
    plot(t*1e3, w_log); grid on;
    legend('w_{iL1}', 'w_{iL2}', 'w_{vCs}', 'w_{vCo}');
    xlabel('t (ms)'); ylabel('\hat{w} (per sample)');
    title('Estimated disturbance');
end

%% ===================== Local functions ==================================
function [A1, B1, C1, A2, B2, C2] = sepic_switched_matrices(p)
% State-space matrices of the two CCM sub-circuits, x' = A x + B [Vin; Vf],
% vo = C x.  1 = switch on / diode off, 2 = switch off / diode on.
    k  = p.R/(p.R + p.rCo);          % output ESR divider
    rx = k*p.rCo + p.Rd;             % resistance seen by (iL1+iL2) in the diode path

    A1 = [-(p.rL1+p.Rds)/p.L1,  -p.Rds/p.L1,                 0,       0;
          -p.Rds/p.L2,          -(p.rCs+p.rL2+p.Rds)/p.L2,   1/p.L2,  0;
           0,                   -1/p.Cs,                     0,       0;
           0,                    0,                          0,      -1/(p.Co*(p.R+p.rCo))];
    B1 = [1/p.L1 0; 0 0; 0 0; 0 0];
    C1 = [0 0 0 k];

    A2 = [-(p.rL1+p.rCs+rx)/p.L1,  -rx/p.L1,            -1/p.L1,  -k/p.L1;
          -rx/p.L2,                -(rx+p.rL2)/p.L2,     0,       -k/p.L2;
           1/p.Cs,                  0,                   0,        0;
           k/p.Co,                  k/p.Co,              0,       -1/(p.Co*(p.R+p.rCo))];
    B2 = [1/p.L1 -1/p.L1; 0 -1/p.L2; 0 0; 0 0];
    C2 = [k*p.rCo, k*p.rCo, 0, k];
end

function [D, x] = sepic_operating_point(A1,B1,C1,A2,B2,C2, Vin, Vf, Vo)
% Duty cycle and averaged steady state that give vo = Vo.
    u   = [Vin; Vf];
    xss = @(d) -(d*A1 + (1-d)*A2) \ ((d*B1 + (1-d)*B2)*u);
    vo  = @(d) (d*C1 + (1-d)*C2) * xss(d);
    % search below the peak of the (non-monotonic) non-ideal gain curve
    dgrid = linspace(0.01, 0.99, 500);
    vgrid = arrayfun(vo, dgrid);
    [vmax, imax] = max(vgrid);
    if vmax < Vo
        error('Vo = %.2f V unreachable at Vin = %.2f V (max %.2f V)', Vo, Vin, vmax);
    end
    D = fzero(@(d) vo(d) - Vo, [dgrid(1) dgrid(imax)]);
    x = xss(D);
end

function x = lti_step(A, B, x, u, h)
% Exact ZOH step of x' = A x + B u over time h.
    n  = size(A,1);
    E  = expm([A B; zeros(size(B,2), n + size(B,2))] * h);
    x  = E(1:n,1:n)*x + E(1:n,n+1:end)*u;
end
