%% Test_X1_wholeleg_v3.m
% Induction Motor Drive Simulation (Scalar V/f + PI speed loop)
% Whole-leg open-circuit fault (phase-A leg open) at tfault.
%
% Key points for WHOLE-LEG fault:
%  - Phase-A current must collapse: ia ≈ 0 after tfault
%  - Remaining phases become two-phase operation: ib ≈ -ic
%  - Torque ripple increases (2*fe component), speed dips then recovers.
%
% This script enforces the physical constraint i_alpha = 0 (equivalent to ia=0
% in abc with floating neutral) by applying the flux-linkage constraint:
%   i_alpha = (Lr*psi_s_alpha - Lm*psi_r_alpha)/Delta = 0
% => psi_s_alpha = (Lm/Lr)*psi_r_alpha
%
% Then it drives ONLY the healthy legs (B and C) via the line-to-line voltage v_bc.
%
% ---------------------------------------------------------------

clear; clc; close all;

%% ================== Machine Parameters (example) ==================
Rs = 0.435;
Rr = 0.816;
Lls = 0.002;          % stator leakage
Llr = 0.002;          % rotor leakage
Lm  = 0.0693;         % magnetizing
Ls  = Lls + Lm;        % total stator
Lr  = Llr + Lm;        % total rotor

P = 4;                % poles
J = 0.089;            % inertia
B = 0.001;            % viscous friction

Delta = Ls*Lr - Lm^2;

%% ================== Simulation Settings ==================
Ts = 1e-4;
Tend = 3.0;
time = 0:Ts:Tend;
N = numel(time);

tfault = 1.0;          % [s] whole-leg open-circuit time

TL = 50;               % [N.m] LOAD TORQUE (your case)

%% ================== Inverter + Scalar Control ==================
Vdc = 600;
Vpk_max = Vdc/sqrt(3); % SVPWM fundamental phase peak limit

% reference speed
w_ref = 1500 * 2*pi/60;     % [rad/s]
ws_ref = (P/2) * w_ref;     % [rad/s] electrical synchronous

% V/f base
f_base   = 50;
ws_base  = 2*pi*f_base;
Vpk_base = 325;             % ~400V L-L rms motor => 325V phase peak
Kvf      = Vpk_base/ws_base;

% Speed PI -> slip compensation (adds to ws_ref)
% (tuned to reduce start-up overshoot vs v2)
Kpw = 0.6;
Kiw = 30.0;

ws_min = 0;
ws_max = 2*pi*80;

%% ================== State Initialization ==================
psi_s = 0 + 1j*0;     % stator flux (alpha + j*beta)
psi_r = 0 + 1j*0;     % rotor  flux (alpha + j*beta)
wm    = 0;            % mech speed [rad/s]

theta = 0;
int_w = 0;

% logs
ia = zeros(1,N); ib = zeros(1,N); ic = zeros(1,N);
wm_log = zeros(1,N);
Te_log = zeros(1,N);

%% ================== Main Loop ==================
for k = 1:N
    t = time(k);

    % ---- electrical rotor speed ----
    omega_r = (P/2)*wm;  % electrical rotor speed [rad/s]

    % ---- WHOLE-LEG FAULT constraint (phase-A open) ----
    fault = (t >= tfault);
    if fault
        % Enforce i_alpha = 0 exactly via flux constraint
        psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
    end

    % ---- currents from flux linkages ----
    i_s = (Lr*psi_s - Lm*psi_r)/Delta;   % stator current (alpha + j beta)
    i_r = (-Lm*psi_s + Ls*psi_r)/Delta;  % rotor current  (alpha + j beta)

    i_alpha = real(i_s);
    i_beta  = imag(i_s);

    % ---- Inverse Clarke -> phase currents (for plotting) ----
    if fault
        % Whole-leg open => ia ≈ 0, ib = -ic (two-phase)
        ia(k) = 0;
        ib(k) =  (sqrt(3)/2)*i_beta;
        ic(k) = -(sqrt(3)/2)*i_beta;
    else
        ia(k) = i_alpha;
        ib(k) = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
        ic(k) = -0.5*i_alpha - (sqrt(3)/2)*i_beta;
    end

    % ---- speed control (slip compensation) ----
    e_w = w_ref - wm;
    int_w = int_w + e_w*Ts;

    wslip = Kpw*e_w + Kiw*int_w;
    ws_unsat = ws_ref + wslip;
    ws_cmd = min(max(ws_unsat, ws_min), ws_max);

    % simple anti-windup
    if ws_cmd ~= ws_unsat
        int_w = int_w + (ws_cmd - ws_unsat)/max(Kiw,1e-9);
    end

    theta = theta + ws_cmd*Ts;

    % ---- V/f voltage magnitude ----
    Vpk = Kvf * ws_cmd;
    Vpk = min(max(Vpk, 0), Vpk_max);

    % ---- commanded balanced phase voltages ----
    va_cmd = Vpk*sin(theta);
    vb_cmd = Vpk*sin(theta - 2*pi/3);
    vc_cmd = Vpk*sin(theta + 2*pi/3);

    % ---- stator voltage in alpha-beta ----
    if ~fault
        % Healthy: balanced 3-phase, floating neutral handled by vn subtraction
        vn = (va_cmd + vb_cmd + vc_cmd)/3;
        va = va_cmd - vn;
        vb = vb_cmd - vn;
        vc = vc_cmd - vn;

        v_alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
        v_beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
        v_s = v_alpha + 1j*v_beta;

    else
        % Whole-leg open (A leg): only B and C legs are driven.
        % Use the line-to-line voltage v_bc to excite the machine.
        % The alpha channel is open (no controlled excitation).
        v_bc = vb_cmd - vc_cmd;                 % line-to-line (b minus c)
        v_alpha = 0;
        v_beta  = (sqrt(3)/3) * v_bc;           % Clarke: v_beta = (√3/3)*(vb - vc)
        v_s = v_alpha + 1j*v_beta;
    end

    % ---- flux dynamics ----
    dpsi_s = v_s - Rs*i_s;
    dpsi_r = -Rr*i_r + 1j*omega_r*psi_r;

    psi_s = psi_s + dpsi_s*Ts;
    psi_r = psi_r + dpsi_r*Ts;

    % Re-apply constraint after integration (keeps i_alpha≈0 tightly)
    if fault
        psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
        i_s = (Lr*psi_s - Lm*psi_r)/Delta;
    end

    % ---- electromagnetic torque ----
    Te = 1.5*(P/2) * imag(conj(psi_s) * i_s);

    % ---- mechanical dynamics ----
    wm = wm + (Te - TL - B*wm)/J * Ts;

    wm_log(k) = wm;
    Te_log(k) = Te;
end

%% ================== Plots ==================
figure; 
plot(time, wm_log*60/(2*pi), 'LineWidth', 1.2); grid on;
xline(tfault,'r--','Fault');
title('Rotor Speed in rpm (Whole-leg open-circuit)'); xlabel('Time [s]'); ylabel('Speed [rpm]');

figure;
subplot(3,1,1); plot(time, ia, 'LineWidth', 1.0); grid on; xline(tfault,'r--','Fault');
title('Stator Phase Currents (Whole-leg open-circuit)'); ylabel('i_a [A]');
subplot(3,1,2); plot(time, ib, 'LineWidth', 1.0); grid on; xline(tfault,'r--','Fault'); ylabel('i_b [A]');
subplot(3,1,3); plot(time, ic, 'LineWidth', 1.0); grid on; xline(tfault,'r--','Fault'); ylabel('i_c [A]'); xlabel('Time [s]');

figure;
plot(time, Te_log, 'LineWidth', 1.2); grid on; xline(tfault,'r--','Fault');
title('Electromagnetic Torque'); xlabel('Time [s]'); ylabel('Torque [N.m]');
