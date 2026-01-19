%% Test_X1_fixed.m
% Induction Motor Drive Simulation with Open-Circuit (Leg) Fault
% Fixes the main issues in Test_X1.m:
%   1) No more DC "pole voltages" with constant switches.
%   2) No more i = v/Rs (this was creating ~690 A DC currents).
%   3) Motor inductances are made physically consistent: Ls = Lls + Lm, Lr = Llr + Lm.
%   4) A simple V/f + PI speed loop is used so speed converges to the reference (instead of ramping).

clear; clc; close all;

%% ================== Simulation Settings ==================
Ts     = 1e-5;      % Time step [s] (no PWM here, but electrical dynamics still require small Ts)
Tend   = 3;         % Total simulation time [s]
tfault = 1;         % Fault injection time [s]
time   = 0:Ts:Tend;
N      = numel(time);

%% ================== Motor Parameters ==================
P  = 4;             % Poles
Rs = 0.435;         % Stator resistance [ohm]
Rr = 0.816;         % Rotor resistance [ohm]

% IMPORTANT:
% In your original file you had: Ls = 0.002, Lr = 0.002, Lm = 0.0693.
% That is physically inconsistent because Lm must be < Ls and < Lr.
% Usually 0.002 H is the LEAKAGE inductance. Total inductance is: Ls = Lls + Lm, Lr = Llr + Lm.
Lls = 0.002;        % Stator leakage inductance [H]
Llr = 0.002;        % Rotor  leakage inductance [H]
Lm  = 0.0693;       % Magnetizing inductance [H]

Ls = Lls + Lm;      % Total stator inductance [H]
Lr = Llr + Lm;      % Total rotor inductance [H]

J  = 0.089;         % Inertia [kg.m^2]
B  = 0.001;         % Viscous friction [N.m.s]
TL = 20;            % Load torque [N.m] (change as needed)

Delta = Ls*Lr - Lm^2;

%% ================== Drive / V-f Settings ==================
Vdc = 600;                      % DC link voltage [V] (for limits only)
Vpk_max = Vdc/2;                % rough phase peak limit
w_ref   = 1500 * 2*pi/60;       % speed reference [rad/s] (1500 rpm)

% Base values for V/f (choose consistent with your intended rating)
f_base  = 50;                   % [Hz] 4-pole synchronous speed ~1500 rpm
w_base  = 2*pi*f_base;          % [rad/s]
Vpk_base = 325;                 % phase peak at base (230 Vrms phase-neutral)

Kvf = Vpk_base / w_base;        % V/f gain

% Simple PI speed controller that adjusts stator frequency (scalar control)
Kpw = 5;
Kiw = 200;
w_min = 0;
w_max = 2*pi*80;                % limit frequency to 80 Hz

%% ================== State Variables (alpha-beta, complex) ==================
psi_s = 0 + 1j*0;               % stator flux (alpha + j beta)
psi_r = 0 + 1j*0;               % rotor  flux (alpha + j beta)
wm    = 0;                      % mechanical speed [rad/s]
theta = 0;                      % stator electrical angle [rad]

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

    % Clarke inverse for phase currents
    ia(k) = i_alpha;
    ib(k) = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
    ic(k) = -0.5*i_alpha - (sqrt(3)/2)*i_beta;

    % ---- Speed PI (frequency command) ----
    e_w   = w_ref - wm;
    int_w = int_w + e_w*Ts;

    ws_cmd = w_base + Kpw*e_w + Kiw*int_w;  % stator electrical frequency [rad/s]
    ws_cmd = min(max(ws_cmd, w_min), w_max);

    theta  = theta + ws_cmd*Ts;

    % ---- V/f voltage magnitude ----
    Vpk = Kvf * ws_cmd;
    Vpk = min(max(Vpk, 0), Vpk_max);

    % ---- ideal phase voltages (line-neutral) ----
    va = Vpk*sin(theta);
    vb = Vpk*sin(theta - 2*pi/3);
    vc = Vpk*sin(theta + 2*pi/3);

    % ---- Open-circuit (LEG) fault approximation ----
    % This emulates an open inverter leg: phase-A is disconnected (va ≈ 0).
    if t >= tfault
        va = 0;
    end

    % Remove zero-sequence to emulate floating neutral: va+vb+vc = 0
    vn = (va + vb + vc)/3;
    va = va - vn;
    vb = vb - vn;
    vc = vc - vn;

    % Clarke transform for stator voltage
    v_alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    v_beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
    v_s = v_alpha + 1j*v_beta;

    % ---- electrical rotor speed ----
    omega_r = (P/2)*wm;   % electrical rotor speed [rad/s]

    % ---- flux dynamics ----
    dpsi_s = v_s - Rs*i_s;
    dpsi_r = -Rr*i_r + 1j*omega_r*psi_r;  % squirrel-cage rotor (v_r = 0)

    psi_s = psi_s + dpsi_s*Ts;
    psi_r = psi_r + dpsi_r*Ts;

    % ---- electromagnetic torque ----
    Te = 1.5*(P/2) * imag(psi_s * conj(i_s));

    % ---- mechanical dynamics ----
    wm = wm + (Te - TL - B*wm)/J * Ts;

    % log
    wm_log(k) = wm;
    Te_log(k) = Te;
end

%% ================== Plots ==================
figure;
plot(time, wm_log, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('Speed [rad/s]');
title('Rotor Speed Response (V/f) with Inverter Open-Circuit (Leg) Fault');
xline(tfault, '--r', 'Fault');

figure;
plot(time, wm_log*60/(2*pi), 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('Speed [rpm]');
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
xlabel('Time [s]');
ylabel('Torque [N.m]');
title('Electromagnetic Torque');
xline(tfault, '--r', 'Fault');
