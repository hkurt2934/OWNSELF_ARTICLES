%% STEP 2 - Healthy Plant: IM in synchronous dq, fed by 2-level VSI (SPWM)
% - You "command" sinusoidal phase voltage references (open-loop)
% - SPWM generates switching states Sa,Sb,Sc
% - Phase voltages are reconstructed ONLY from (Sa,Sb,Sc,Vdc)
% - Then abc->dq and the IM dq plant is integrated
%
% NOTE: Still NO control (no FOC, no MRAS). This is healthy inverter-fed baseline.

clear; clc; close all;

%% -------------------- Motor parameters (typical 1.5 kW, 4-pole) ----------
p   = 2;            % pole pairs
Rs  = 4.10;         % [ohm]
Rr  = 3.90;         % [ohm]
Lls = 8.0e-3;       % [H]
Llr = 8.0e-3;       % [H]
Lm  = 0.135;        % [H]

J   = 0.008;        % [kg.m^2]
B   = 5e-4;         % [N.m.s]

Ls  = Lls + Lm;
Lr  = Llr + Lm;

% Inductance matrix for dq currents -> dq fluxes (power-invariant form)
Lmat = [ Ls   0   Lm   0;
          0  Ls    0  Lm;
         Lm   0   Lr   0;
          0  Lm    0  Lr ];

%% -------------------- Electrical fundamental (50 Hz) ---------------------
f1  = 50;                 % fundamental frequency [Hz]
we  = 2*pi*f1;             % synchronous electrical speed [rad/s]

%% -------------------- DC link voltage (choose one) -----------------------
Vdc = 560;                 % [V] typical from 400V rectified (~540-580V)
% Vdc = 600;               % [V] also common in simulations

%% -------------------- Inverter modulation (SPWM) -------------------------
fsw = 10e3;                % switching frequency [Hz]
% Modulation index must be <= 1 for linear SPWM.
% The fundamental phase voltage peak achievable (approx) is ~ (m*Vdc/2).
m   = 0.95;                % modulation index (0..1), keep <1

% Fundamental phase voltage reference (peak) for SPWM
Vph_ref_pk = m*(Vdc/2);    % [V] (SPWM linear region)

%% -------------------- Load torque ---------------------------------------
P_rated = 1.5e3;
n_sync  = 1500;
w_sync  = 2*pi*n_sync/60;
T_rated = P_rated / w_sync;   % ~9.55 N.m

Tl = T_rated;   % baseline rated load

%% -------------------- Simulation ----------------------------------------
Ts   = 1e-4;     % 100 us (same as your plan)
Tend = 3.0;      % [s]
t = (0:Ts:Tend).';
N = numel(t);

% States: x = [ids iqs idr iqr wm]^T
x = zeros(5,1);

% Logs
rpm = zeros(N,1);
Te  = zeros(N,1);
ids = zeros(N,1); iqs = zeros(N,1);
va = zeros(N,1); vb = zeros(N,1); vc = zeros(N,1);
Sa_log = zeros(N,1); Sb_log = zeros(N,1); Sc_log = zeros(N,1);

%% -------------------- Main loop -----------------------------------------
for k = 1:N
    tk = t(k);

    % Electrical angle for fundamental
    theta1 = we*tk;

    % --- Open-loop sinusoidal references (what your modulator wants) ---
    % These are NOT applied directly to the motor.
    va_ref = Vph_ref_pk*sin(theta1);
    vb_ref = Vph_ref_pk*sin(theta1 - 2*pi/3);
    vc_ref = Vph_ref_pk*sin(theta1 + 2*pi/3);

    % --- Convert references to duty ratios (centered) ---
    % For a 2-level inverter with pole voltage +/-Vdc/2, a simple mapping is:
    % v_ref = (2*d - 1)*(Vdc/2) => d = 0.5 + v_ref/Vdc
    da = 0.5 + va_ref/Vdc;
    db = 0.5 + vb_ref/Vdc;
    dc = 0.5 + vc_ref/Vdc;

    % Clamp duties to [0,1]
    da = min(max(da,0),1);
    db = min(max(db,0),1);
    dc = min(max(dc,0),1);

    % --- Reconstruct inverter output phase voltages from (S,Vdc) ---
    
    % Convert to phase voltages w.r.t. floating neutral (remove common-mode):
    % --- Average inverter model (NO switching) ---
    va_k = (2*da - 1)*(Vdc/2);
    vb_k = (2*db - 1)*(Vdc/2);
    vc_k = (2*dc - 1)*(Vdc/2);
    
    % Remove common-mode (floating neutral)
    v_cm = (va_k + vb_k + vc_k)/3;
    va_k = va_k - v_cm;
    vb_k = vb_k - v_cm;
    vc_k = vc_k - v_cm;

    % --- Synchronous dq frame angle for plant equations ---
    % Here we use the same synchronous frame as Step-1: theta_e = we*t
    theta_e = theta1;

    % abc -> dq for applied stator voltages
    [vds, vqs] = abc_to_dq(va_k, vb_k, vc_k, theta_e);

    % --- Unpack states ---
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_k  = x(5);

    % Rotor electrical speed and slip
    wr_e = p*wm_k;
    wsl  = we - wr_e;

    % Fluxes from currents
    i4  = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat*i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);

    % Torque
    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k);

    % Flux derivatives in synchronous frame
    dpsi = zeros(4,1);
    dpsi(1) = vds - Rs*ids_k + we*psiqs;
    dpsi(2) = vqs - Rs*iqs_k - we*psids;
    dpsi(3) =      - Rr*idr_k + wsl*psiqr;
    dpsi(4) =      - Rr*iqr_k - wsl*psidr;

    % Convert dpsi -> di
    di = Lmat \ dpsi;

    % Mechanical
    dwm = (Te_k - Tl - B*wm_k)/J;

    % Integrate
    x(1:4) = x(1:4) + Ts*di;
    x(5)   = x(5)   + Ts*dwm;

    % Logs
    rpm(k) = wm_k*60/(2*pi);
    Te(k)  = Te_k;
    ids(k) = ids_k; iqs(k) = iqs_k;

    va(k) = va_k; vb(k) = vb_k; vc(k) = vc_k;
end

%% -------------------- Plots ---------------------------------------------
figure; plot(t, rpm, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title(sprintf('Healthy IM (VSI-SPWM, Vdc=%.0fV, fsw=%.0fkHz), Load=%.2f N·m', Vdc, fsw/1e3, Tl));

figure; plot(t, Te, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Torque [N·m]');
title('Electromagnetic Torque');

figure; plot(t, ids, 'LineWidth', 1.2); hold on;
plot(t, iqs, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Current [A]');
title('Stator Currents in synchronous dq');
legend('i_{ds}','i_{qs}');

% Optional: show a short window of phase voltage PWM (zoom first 5 ms)
figure;
idx = t <= 0.005;
plot(t(idx), va(idx), 'LineWidth', 1.0); hold on;
plot(t(idx), vb(idx), 'LineWidth', 1.0);
plot(t(idx), vc(idx), 'LineWidth', 1.0);
grid on; xlabel('Time [s]'); ylabel('Voltage [V]');
title('Inverter phase voltages (PWM waveform) - first 5 ms');
legend('v_a','v_b','v_c');

%% -------------------- Local functions -----------------------------------

function [vd, vq] = abc_to_dq(va, vb, vc, theta)
    % Power-invariant Clarke
    alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
    % Park
    c = cos(theta); s = sin(theta);
    vd =  alpha*c + beta*s;
    vq = -alpha*s + beta*c;
end
