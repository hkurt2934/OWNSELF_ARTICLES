%% STEP 2 (Healthy) - IM Plant in synchronous dq, fed by 2-level VSI using SVPWM (average model)
% - No control (no FOC, no MRAS)
% - 50 Hz open-loop voltage vector reference
% - SVPWM creates duty ratios da,db,dc
% - Average phase voltages are reconstructed using only Vdc and duties
% - Plant is induction motor in synchronous dq frame (currents as states)

clear; clc; close all;

%% -------------------- Motor parameters (typical 1.5 kW, 4-pole) ----------
p   = 2;            % pole pairs
Rs  = 4.10;         % stator resistance [ohm]
Rr  = 3.90;         % rotor resistance [ohm]
Lls = 8.0e-3;       % stator leakage inductance [H]
Llr = 8.0e-3;       % rotor leakage inductance [H]
Lm  = 0.135;        % magnetizing inductance [H]

J   = 0.008;        % inertia [kg.m^2]
B   = 5e-4;         % viscous friction [N.m.s]

Ls  = Lls + Lm;
Lr  = Llr + Lm;

% Inductance matrix for dq currents -> dq fluxes
% i = [ids iqs idr iqr]^T, psi = Lmat*i
Lmat = [ Ls   0   Lm   0;
          0  Ls    0  Lm;
         Lm   0   Lr   0;
          0  Lm    0  Lr ];

%% -------------------- Electrical fundamental (50 Hz) ---------------------
f1  = 50;                 % [Hz]
we  = 2*pi*f1;             % synchronous electrical speed [rad/s]
theta0 = 0;

%% -------------------- DC link voltage -----------------------------------
Vdc = 560;                 % [V] typical from 400V rectified (540-580V)
% Vdc = 600;               % [V] also fine

%% -------------------- Open-loop voltage level (SVPWM linear region) ------
% Desired motor rating: ~400 V line-line RMS, 50 Hz.
% Required phase peak for 400V_LL,rms:
Vll_rms_des = 400;
Vph_pk_req  = sqrt(2)*(Vll_rms_des/sqrt(3));   % ~326.6 V

% SVPWM linear max phase peak is about Vdc/sqrt(3) (~0.577*Vdc)
Vph_pk_max  = Vdc/sqrt(3);

% Choose reference safely inside SVPWM linear limit
Vph_pk_ref  = min(Vph_pk_req, 0.98*Vph_pk_max);

fprintf('Requested Vph_pk for 400VLL: %.1f V\n', Vph_pk_req);
fprintf('SVPWM max Vph_pk (linear):   %.1f V\n', Vph_pk_max);
fprintf('Using Vph_pk_ref:            %.1f V\n', Vph_pk_ref);

%% -------------------- Load torque ---------------------------------------
P_rated = 1.5e3;
n_sync  = 1500;
w_sync  = 2*pi*n_sync/60;
T_rated = P_rated / w_sync;      % ~9.55 N.m
Tl = T_rated;                    % rated load baseline

%% -------------------- Simulation settings --------------------------------
Ts   = 1e-4;       % 100 us
Tend = 3.0;        % [s]
t = (0:Ts:Tend).';
N = numel(t);

% States: x = [ids iqs idr iqr wm]^T
x = zeros(5,1);

% Logs
rpm = zeros(N,1);
Te  = zeros(N,1);
ids = zeros(N,1); iqs = zeros(N,1);
idr = zeros(N,1); iqr = zeros(N,1);
va  = zeros(N,1); vb  = zeros(N,1); vc  = zeros(N,1);
vds_log = zeros(N,1); vqs_log = zeros(N,1);
da_log = zeros(N,1); db_log = zeros(N,1); dc_log = zeros(N,1);

%% -------------------- Main loop -----------------------------------------
for k = 1:N
    tk = t(k);

    % Synchronous electrical angle used by dq plant and open-loop reference
    theta_e = theta0 + we*tk;

    % --- Open-loop voltage vector reference in alpha-beta ---
    % If you want a pure rotating vector at 50 Hz:
    v_alpha_ref = Vph_pk_ref*cos(theta_e);
    v_beta_ref  = Vph_pk_ref*sin(theta_e);

    % --- SVPWM duties from v_alpha_ref, v_beta_ref, Vdc ---
    [da, db, dc] = svpwm_duty(v_alpha_ref, v_beta_ref, Vdc);

    % --- Average phase voltages applied to motor (floating neutral) ---
    [va_k, vb_k, vc_k] = duty_to_vabc(da, db, dc, Vdc);

    % abc -> dq (applied stator voltages)
    [vds, vqs] = abc_to_dq(va_k, vb_k, vc_k, theta_e);

    % Unpack states
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_k  = x(5);

    % Rotor electrical speed and slip frequency
    wr_e = p*wm_k;
    wsl  = we - wr_e;

    % Fluxes
    i4  = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat*i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);

    % Electromagnetic torque
    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k);

    % dq flux derivatives in synchronous frame
    dpsi = zeros(4,1);
    dpsi(1) = vds - Rs*ids_k + we*psiqs;
    dpsi(2) = vqs - Rs*iqs_k - we*psids;
    dpsi(3) =      - Rr*idr_k + wsl*psiqr;
    dpsi(4) =      - Rr*iqr_k - wsl*psidr;

    % Convert dpsi -> di  (Lmat * di = dpsi)
    di = Lmat \ dpsi;

    % Mechanical dynamics
    dwm = (Te_k - Tl - B*wm_k)/J;

    % Integrate (Euler)
    x(1:4) = x(1:4) + Ts*di;
    x(5)   = x(5)   + Ts*dwm;

    % Logs
    rpm(k) = wm_k*60/(2*pi);
    Te(k)  = Te_k;
    ids(k) = ids_k; iqs(k) = iqs_k; idr(k)=idr_k; iqr(k)=iqr_k;

    va(k) = va_k; vb(k) = vb_k; vc(k) = vc_k;
    vds_log(k)=vds; vqs_log(k)=vqs;
    da_log(k)=da; db_log(k)=db; dc_log(k)=dc;
end

%% -------------------- Plots ---------------------------------------------
figure; plot(t, rpm, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title(sprintf('Healthy IM (VSI-SVPWM avg, Vdc=%.0f V), Load=%.2f N·m', Vdc, Tl));

figure; plot(t, Te, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Torque [N·m]');
title('Electromagnetic Torque');

figure; plot(t, ids, 'LineWidth', 1.2); hold on;
plot(t, iqs, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Current [A]');
title('Stator Currents in synchronous dq');
legend('i_{ds}','i_{qs}','Location','best');

% Optional: show duties for sanity (first 40 ms)
figure;
idx = t <= 0.04;
plot(t(idx), da_log(idx), 'LineWidth', 1.1); hold on;
plot(t(idx), db_log(idx), 'LineWidth', 1.1);
plot(t(idx), dc_log(idx), 'LineWidth', 1.1);
grid on; xlabel('Time [s]'); ylabel('Duty');
title('SVPWM duty ratios (first 40 ms)');
legend('d_a','d_b','d_c','Location','best');

% Optional: show applied phase voltages (first 5 ms)
figure;
idx2 = t <= 0.005;
plot(t(idx2), va(idx2), 'LineWidth', 1.0); hold on;
plot(t(idx2), vb(idx2), 'LineWidth', 1.0);
plot(t(idx2), vc(idx2), 'LineWidth', 1.0);
grid on; xlabel('Time [s]'); ylabel('Voltage [V]');
title('Average phase voltages v_a,v_b,v_c (first 5 ms)');
legend('v_a','v_b','v_c','Location','best');

%% -------------------- Local functions -----------------------------------
function [da, db, dc] = svpwm_duty(v_alpha, v_beta, Vdc)
% SVPWM duty computation (average model, 2-level VSI).
% Input v_alpha, v_beta are desired phase-to-neutral equivalent alpha-beta voltages [V].
% Outputs da,db,dc in [0,1].

    % Linear SVPWM limit for alpha-beta vector magnitude
    Vmax = Vdc/sqrt(3);
    Vref = hypot(v_alpha, v_beta);
    if Vref > Vmax
        scale = Vmax / Vref;
        v_alpha = v_alpha * scale;
        v_beta  = v_beta  * scale;
    end

    % Compute equivalent three-phase "commands" (normalized) from alpha-beta
    % (power-invariant mapping consistent with Clarke used elsewhere)
    ua = (2/3)*( v_alpha )/Vdc;
    ub = (2/3)*(-v_alpha/2 + (sqrt(3)/2)*v_beta)/Vdc;
    uc = (2/3)*(-v_alpha/2 - (sqrt(3)/2)*v_beta)/Vdc;

    % Zero-sequence injection (equivalent SVPWM)
    u_max = max([ua, ub, uc]);
    u_min = min([ua, ub, uc]);
    u0 = -0.5*(u_max + u_min);

    da = 0.5 + (ua + u0);
    db = 0.5 + (ub + u0);
    dc = 0.5 + (uc + u0);

    da = min(max(da,0),1);
    db = min(max(db,0),1);
    dc = min(max(dc,0),1);
end

function [va, vb, vc] = duty_to_vabc(da, db, dc, Vdc)
% Duty ratios -> average phase voltages (floating neutral)
    v_aO = (2*da - 1)*(Vdc/2);
    v_bO = (2*db - 1)*(Vdc/2);
    v_cO = (2*dc - 1)*(Vdc/2);

    v_cm = (v_aO + v_bO + v_cO)/3;
    va = v_aO - v_cm;
    vb = v_bO - v_cm;
    vc = v_cO - v_cm;
end

function [vd, vq] = abc_to_dq(va, vb, vc, theta)
% abc -> dq using power-invariant Clarke + Park
    alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
    c = cos(theta); s = sin(theta);
    vd =  alpha*c + beta*s;
    vq = -alpha*s + beta*c;
end
