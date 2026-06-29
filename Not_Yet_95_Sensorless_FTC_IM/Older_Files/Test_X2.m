%% Test_X1_fixed_v2.m
% Induction Motor Drive Simulation with Inverter Open-Circuit (Leg) Fault
%
% What this v2 fixes (based on your latest plots):
%  1) Torque sign was inverted -> speed ran NEGATIVE. (Main root cause)
%  2) Stator frequency command referenced incorrectly -> unnecessary saturation.
%  3) Voltage limit updated for SVPWM (V_phase,peak,max = Vdc/sqrt(3)).
%
% Notes:
%  - This is still a *scalar* V/f + PI speed loop (not FOC). It is meant to
%    generate physically plausible currents/torque/speed and show the effect
%    of a leg open-circuit fault.
%  - A true IGBT open-circuit (single device) needs diode-current-direction
%    logic or a switching model (Simulink Universal Bridge). If you want
%    that, tell me which device (Sa+, Sa-, Sb+, ...).

clear; clc; close all;

%% ================== Simulation Settings ==================
Ts     = 1e-5;      % Time step [s]
Tend   = 3;         % Total simulation time [s]
tfault = 1;         % Fault injection time [s]
time   = 0:Ts:Tend;
N      = numel(time);

%% ================== Motor Parameters ==================
P  = 4;             % Poles
Rs = 0.435;         % Stator resistance [ohm]
Rr = 0.816;         % Rotor resistance [ohm]

% Inductances (consistent): Ls = Lls + Lm, Lr = Llr + Lm
Lls = 0.002;        % Stator leakage [H]
Llr = 0.002;        % Rotor  leakage [H]
Lm  = 0.0693;       % Magnetizing [H]

Ls = Lls + Lm;
Lr = Llr + Lm;

J  = 0.089;         % Inertia [kg.m^2]
B  = 0.001;         % Viscous friction [N.m.s]
TL = 20;            % Load torque [N.m] (increase/decrease to match your test)

Delta = Ls*Lr - Lm^2;

%% ================== Inverter / V-f Settings ==================
Vdc = 600;                          % DC link [V]
Vpk_max = Vdc/sqrt(3);              % SVPWM fundamental phase peak limit ~0.577*Vdc

% Mechanical speed reference
w_ref = 1500 * 2*pi/60;             % [rad/s] 1500 rpm

% Electrical stator speed reference (approx synchronous)
ws_ref = (P/2) * w_ref;             % [rad/s]

% Base V/f (choose consistent with your intended rating)
f_base   = 50;                      % [Hz]
ws_base  = 2*pi*f_base;             % [rad/s]
Vpk_base = 325;                     % [V] phase peak at base
Kvf      = Vpk_base/ws_base;        % [V/(rad/s)]

% PI speed controller output = slip compensation (added to ws_ref)
Kpw = 1.0;
Kiw = 50.0;

ws_min = 0;
ws_max = 2*pi*80;                   % [rad/s]

%% ================== States (alpha-beta, complex) ==================
psi_s = 0 + 1j*0;                   % stator flux
psi_r = 0 + 1j*0;                   % rotor  flux
wm    = 0;                          % mechanical speed [rad/s]
theta = 0;                          % stator electrical angle [rad]

int_w = 0;

%% ================== Logging ==================
wm_log = zeros(1,N);
Te_log = zeros(1,N);
ia = zeros(1,N); ib = zeros(1,N); ic = zeros(1,N);

%% ================== Main Loop ==================
for k = 1:N
    t = time(k);

    % ---- currents from flux linkages ----
    i_s = (Lr*psi_s - Lm*psi_r)/Delta;      % stator current (alpha + j beta)
    i_r = (-Lm*psi_s + Ls*psi_r)/Delta;     % rotor current  (alpha + j beta)

    i_alpha = real(i_s);
    i_beta  = imag(i_s);

    % Inverse Clarke -> phase currents
    ia(k) = i_alpha;
    ib(k) = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
    ic(k) = -0.5*i_alpha - (sqrt(3)/2)*i_beta;

    % ---- Speed PI (adds slip compensation to ws_ref) ----
    e_w = w_ref - wm;

    % simple anti-windup: integrate only if not saturating hard
    int_w = int_w + e_w*Ts;
    ws_unsat = ws_ref + Kpw*e_w + Kiw*int_w;
    ws_cmd   = min(max(ws_unsat, ws_min), ws_max);
    if ws_cmd ~= ws_unsat
        % back-calculate one-step correction (keeps integrator bounded)
        int_w = int_w + (ws_cmd - ws_unsat)/max(Kiw,1e-9);
    end

    theta = theta + ws_cmd*Ts;

    % ---- V/f voltage magnitude ----
    Vpk = Kvf * ws_cmd;
    Vpk = min(max(Vpk, 0), Vpk_max);

    % ---- ideal phase voltages (line-neutral) ----
    va = Vpk*sin(theta);
    vb = Vpk*sin(theta - 2*pi/3);
    vc = Vpk*sin(theta + 2*pi/3);

    % ---- Open-circuit (LEG) fault approximation ----
    % Leg open-circuit: phase-A terminal is disconnected from DC bus.
    % A quick approximation is forcing the phase-A inverter voltage to 0.
    if t >= tfault
        va = 0;
    end

    % Floating neutral (remove zero-sequence): va+vb+vc = 0
    vn = (va + vb + vc)/3;
    va = va - vn;
    vb = vb - vn;
    vc = vc - vn;

    % Clarke -> alpha-beta stator voltage
    v_alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    v_beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
    v_s = v_alpha + 1j*v_beta;

    % ---- electrical rotor speed ----
    omega_r = (P/2)*wm;

    % ---- flux dynamics ----
    dpsi_s = v_s - Rs*i_s;
    dpsi_r = -Rr*i_r + 1j*omega_r*psi_r;   % stationary frame rotor eq.

    psi_s = psi_s + dpsi_s*Ts;
    psi_r = psi_r + dpsi_r*Ts;

    % ---- electromagnetic torque (SIGN FIX) ----
    % Using T = (3/2)*(P/2)*(psi_alpha*i_beta - psi_beta*i_alpha)
    % which equals imag(conj(psi_s)*i_s)
    Te = 1.5*(P/2) * imag(conj(psi_s) * i_s);

    % ---- mechanical dynamics ----
    wm = wm + (Te - TL - B*wm)/J * Ts;

    wm_log(k) = wm;
    Te_log(k) = Te;
end

%% ================== Plots ==================
figure;
plot(time, wm_log, 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Speed [rad/s]');
title('Rotor Speed Response (V/f) with Inverter Open-Circuit (Leg) Fault');
xline(tfault, '--r', 'Fault');

figure;
plot(time, wm_log*60/(2*pi), 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title('Rotor Speed in rpm');
xline(tfault, '--r', 'Fault');

figure;
subplot(3,1,1);
plot(time, ia, 'LineWidth', 1.0); 
hold on;
legend('i_a','Location','best');
xlabel('Time [s]');
ylabel('Phase Currents [A]');
title('Stator Phase Currents under Inverter Open-Circuit (Leg) Fault');
subplot(3,1,2);
plot(time, ib, 'LineWidth', 1.0);
hold on;
legend('i_b','Location','best');
xlabel('Time [s]');
ylabel('Phase Currents [A]');
title('Stator Phase Currents under Inverter Open-Circuit (Leg) Fault');
subplot(3,1,3);
plot(time, ic, 'LineWidth', 1.0);
hold on;
legend('i_c','Location','best');
xlabel('Time [s]');
ylabel('Phase Currents [A]');
title('Stator Phase Currents under Inverter Open-Circuit (Leg) Fault');
xline(tfault, '--r', 'Fault');
grid on;

figure;
plot(time, Te_log, 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Torque [N.m]');
title('Electromagnetic Torque');
xline(tfault, '--r', 'Fault');
