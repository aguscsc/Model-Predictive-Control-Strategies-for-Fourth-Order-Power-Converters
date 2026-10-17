%% H-infinity control for the NON-IDEAL SEPIC
% Same plant, parasitics, reference generator and simulation as
% sepic_mpc_nonideal.m, but the controller is a discrete-time H-infinity
% state-feedback law with integral action on vo:
%
%     d(k) = D_ref + K * [x(k) - x_ref ; xi(k)],   xi(k+1) = xi(k) + Ts*(vo - Vo)
%
% K is synthesised by solving the discrete game Riccati equation (idare
% with an indefinite weight) and bisecting on gamma, so only the Control
% System Toolbox is needed (no Robust Control Toolbox / hinfsyn).
%
% Generalized plant used for the synthesis (deviation coordinates):
%   xa(k+1) = Aa xa + Bu du + Bw w        xa = [dx; xi]
%   z       = [Wx^(1/2) xa ; wu du]       performance output
%   w       = [dVin ; dIo]                input-voltage and load-current disturbances
% The controller minimises gamma such that ||T_zw||_inf < gamma.
%
% States:  x = [iL1; iL2; vCs; vCo]  (vCs, vCo behind the ESR; vo is an output)
% CCM is assumed throughout.
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

% --- H-inf performance weights (Bryson-style tolerances) ---
% Tuned for the smoothest vo: after a Vin or load step vo makes a single
% hump and recovers like a first-order system (tau ~ 1 ms), no overshoot.
% To get that, the current and vCs weights are made very loose (10 A,
% 10 A, 30 V tolerances = practically unweighted), so the controller only
% acts on vo, its integral and d.
% TRADE-OFF: the Cs/L1/L2 resonance (~2.5 kHz) is then left almost
% undamped (zeta ~0.02): iL1, iL2 and vCs ring for >10 ms after a step,
% even though vo does not. Previous tuning (damps that ring, but vo
% oscillates) in brackets.
diL1_max  = 10;     % fraction of nominal iL1                     [0.037]
diL2_max  = 10;     % fraction of nominal iL2                     [0.24]
dvCs_max  = 30;     % V                                           [0.22]
dvo_max   = 0.05;   % V  (kept fixed: sets the overall scale)
dmode_max = 1;      % A, tolerance on the internal-mode imbalance [0.017]
w_mode    = 1;      % 1 = internal-mode term on, 0 = off
T_int     = 1e-3;   % s, integral time scale: xi tolerance = dvo_max*T_int [0.22e-3]
dd_max    = 0.03;   % duty-cycle deviation tolerance              [0.0027]

% --- H-inf disturbance scaling (expected size of each disturbance) ---
dVin_max  = 1.0;    % V, input-voltage disturbance (kept fixed: sets the scale)
dIo_max   = 0.2;    % A, load-current disturbance injected at the output [0.45]

% --- gamma bisection ---
gamma_lo  = 1e-4;   % lower bound for the search
gamma_hi  = 1e4;    % upper bound for the search
gamma_tol = 1e-3;   % relative tolerance
gamma_backoff = 1.01;% use gamma = backoff*gamma_opt (central controller near
                    % gamma_opt has very high gain; >1 gives a softer K) [1.03]

% --- Constraints / anti-windup ---
dmin = 0.05;
dmax = 0.75;

% --- Reference generator ---
alpha      = 0.01;  % Vin low-pass filter coefficient
Vin_lut    = 8:0.05:16;  % Vin grid for the steady-state reference table

% --- Simulation ---
Tsim         = 8e-3;               % s
t_vin_step   = 2e-3;                % s, time of the Vin step
Vin_step     = 14;                % V, Vin after the step
t_load_step  = 5e-3;              % s, time of the load step (Inf = no step)
R_step       = 30;                  % ohm, load after the step
x0_kick      = [0; -0; 0; 0];   % initial deviation from the operating point
use_ctrl     = true;                % false = open loop (d held at nominal)
plant_model  = 'averaged';          % 'averaged' or 'switched' (on/off sub-intervals)
    
%% Pack parameters
par = struct('L1',L1,'L2',L2,'Cs',Cs,'Co',Co,'R',R, ...
             'rL1',rL1+Rsh1,'rL2',rL2+Rsh2,'rCs',rCs,'rCo',rCo, ...
             'Rds',Rds_on,'Vf',Vf,'Rd',Rd);
[A1, B1, C1, A2, B2, C2] = sepic_switched_matrices(par);

%% Non-ideal operating point: find D that gives vo = Vo at Vin
[D_nom, x_nom] = sepic_operating_point(A1,B1,C1,A2,B2,C2, Vin, Vf, Vo);
il1 = x_nom(1);  il2 = x_nom(2);  vcs = x_nom(3);  vco = x_nom(4);

Pin  = Vin*il1;
Pout = Vo^2/R;
fprintf('Operating point (non-ideal):\n');
fprintf('  D = %.4f (ideal %.4f)\n', D_nom, Vo/(Vo + Vin));
fprintf('  iL1 = %.4f A, iL2 = %.4f A, vCs = %.3f V, vCo = %.3f V\n', il1, il2, vcs, vco);
fprintf('  efficiency = %.1f %%\n', 100*Pout/Pin);

%% Small-signal model around (D_nom, x_nom), ZOH discretised
U_nom = [Vin; Vf];
A_hat = D_nom*A1 + (1-D_nom)*A2;
Bd_c  = (A1 - A2)*x_nom + (B1 - B2)*U_nom;          % d   -> x
Bv_c  = D_nom*B1(:,1) + (1-D_nom)*B2(:,1);          % Vin -> x
Bi_c  = [0; 0; 0; -1/Co];                           % io  -> x (extra load current)

sys_x = ss(A_hat, [Bd_c Bv_c Bi_c], eye(4), zeros(4,3));
sysd  = c2d(sys_x, Ts, 'zoh');
Ad = sysd.A;
Bu = sysd.B(:,1);
Bw = sysd.B(:,2:3) * diag([dVin_max, dIo_max]);     % scaled disturbances

%% Augment with an integrator on the vo error (vCo state)
nx = 4;
Aa  = [Ad zeros(nx,1); Ts*[0 0 0 1] 1];
Bua = [Bu; 0];
Bwa = [Bw; zeros(1,2)];
na  = nx + 1;
nw  = size(Bwa, 2);

%% Performance weights
c_mode = [1 -1 0 0];
Qx     = diag(1 ./ [diL1_max*il1, diL2_max*il2, dvCs_max, dvo_max].^2) ...
       + w_mode * (c_mode' * c_mode) / dmode_max^2;
Qa     = blkdiag(Qx, 1/(dvo_max*T_int)^2);
Ru     = 1 / dd_max^2;

%% gamma bisection on the discrete H-inf state-feedback Riccati equation
if ~hinf_sf_feasible(Aa, Bua, Bwa, Qa, Ru, gamma_hi)
    error('H-inf problem infeasible at gamma_hi = %g; raise gamma_hi or relax the weights', gamma_hi);
end
lo = gamma_lo;  hi = gamma_hi;
while hi/lo > 1 + gamma_tol
    g = sqrt(lo*hi);                                 % bisect in log scale
    if hinf_sf_feasible(Aa, Bua, Bwa, Qa, Ru, g)
        hi = g;
    else
        lo = g;
    end
end
gamma_opt = hi;
gamma     = gamma_backoff * gamma_opt;
[~, K]    = hinf_sf_feasible(Aa, Bua, Bwa, Qa, Ru, gamma);   % du = -K*xa

%% Verify the design
Acl  = Aa - Bua*K;
Cz   = [sqrtm(Qa); -sqrt(Ru)*K];
Tzw  = ss(Acl, Bwa, Cz, zeros(size(Cz,1), nw), Ts);
gamma_achieved = norm(Tzw, inf);
fprintf('gamma_opt = %.4g | design gamma = %.4g | achieved ||Tzw||inf = %.4g\n', ...
        gamma_opt, gamma, gamma_achieved);
fprintf('K = [%s]\n', sprintf(' %.4g', K));
pcl = eig(Acl);
fprintf('max |closed-loop pole| = %.6f\n', max(abs(pcl)));
disp('Closed-loop poles (continuous equivalent, rad/s):');
disp(log(pcl)/Ts);

% Closed-loop disturbance -> vo responses
Cvo = [0 0 0 1 0];
T_vin = ss(Acl, Bwa(:,1)/dVin_max, Cvo, 0, Ts);
T_io  = ss(Acl, Bwa(:,2)/dIo_max,  Cvo, 0, Ts);
figure;
bodemag(T_vin, T_io); grid on;
legend('v_{in} \rightarrow v_o', 'i_o \rightarrow v_o (ohm)');
title('Closed-loop disturbance rejection (H_\infty)');

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

x        = x_nom + x0_kick;      % physical (large-signal) plant state
xi       = 0;                    % integrator state
vin_filt = Vin;
vo_meas  = (D_nom*C1 + (1-D_nom)*C2) * x;

x_log    = zeros(4, N);
xref_log = zeros(4, N);
vo_log   = zeros(1, N);
d_log    = zeros(1, N);

for k = 1:N
    % --- measure ---
    vin_meas = vin_profile(k);

    % --- adaptive reference from the non-ideal steady-state table ---
    vin_filt = alpha*vin_meas + (1-alpha)*vin_filt;
    D_ref_k  = interp1(Vin_lut, D_lut, vin_filt, 'linear', 'extrap');
    x_ref_k  = interp1(Vin_lut, X_lut.', vin_filt, 'linear', 'extrap').';

    % --- controller ---
    if use_ctrl
        d_unsat = D_ref_k - K * [x - x_ref_k; xi];
        d_k     = min(max(d_unsat, dmin), dmax);
        if d_k == d_unsat                          % anti-windup: freeze xi when saturated
            xi = xi + Ts*(vo_meas - Vo);
        end
    else
        d_k = D_nom;
    end

    % --- plant (load may change; rebuild sub-circuit matrices if it did) ---
    if k == 1 || R_profile(k) ~= R_profile(k-1)
        par.R = R_profile(k);
        [A1p, B1p, C1p, A2p, B2p, C2p] = sepic_switched_matrices(par);
    end

    % --- log (vo uses the averaged output equation) ---
    vo_meas       = (d_k*C1p + (1-d_k)*C2p) * x;
    x_log(:,k)    = x;
    xref_log(:,k) = x_ref_k;
    vo_log(k)     = vo_meas;
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
names = {'i_{L1} (A)', 'i_{L2} (A)', 'v_{Cs} (V)', 'v_o (V)'};
figure;
for n = 1:3
    subplot(5,1,n);
    plot(t*1e3, x_log(n,:), t*1e3, xref_log(n,:), '--'); grid on;
    ylabel(names{n});
end
subplot(5,1,4);
plot(t*1e3, vo_log, t*1e3, Vo*ones(1,N), '--'); grid on;
ylabel(names{4});
subplot(5,1,5);
plot(t*1e3, d_log); grid on;
ylabel('d'); xlabel('t (ms)');
sgtitle(sprintf('Non-ideal SEPIC H_\\infty (%s plant)', plant_model));

%% ===================== Local functions ==================================
function [ok, K] = hinf_sf_feasible(A, Bu, Bw, Q, R, gamma)
% Discrete-time H-inf state feedback for z = [Q^(1/2) x; R^(1/2) u]:
% solve the game Riccati equation
%   X = A'XA + Q - A'X[Bu Bw] (blkdiag(R,-gamma^2 I) + [Bu Bw]'X[Bu Bw])^-1 [Bu Bw]'XA
% and check X >= 0, gamma^2 I - Bw'XBw > 0 and closed-loop stability.
% Returns the central state-feedback gain (u = -K x) when feasible.
    K  = [];
    nw = size(Bw, 2);
    Bt = [Bu Bw];
    Rt = blkdiag(R, -gamma^2*eye(nw));
    try
        [X, ~, ~, info] = idare(A, Bt, Q, Rt);
    catch
        ok = false;  return;
    end
    if isempty(X) || info.Report > 1 || any(~isfinite(X(:)))
        ok = false;  return;
    end
    X = (X + X')/2;
    if min(eig(X)) < -1e-9*max(1, norm(X))
        ok = false;  return;
    end
    G = gamma^2*eye(nw) - Bw'*X*Bw;
    if min(eig((G + G')/2)) <= 0
        ok = false;  return;
    end
    Xb = X + X*Bw*(G \ (Bw'*X));                    % worst-case-disturbance-adjusted cost
    K  = (R + Bu'*Xb*Bu) \ (Bu'*Xb*A);
    ok = max(abs(eig(A - Bu*K))) < 1;
end

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
