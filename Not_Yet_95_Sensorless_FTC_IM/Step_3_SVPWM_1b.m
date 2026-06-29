
%% STEP 3 (Healthy) - Induction Motor FOC (IFOC) + SVPWM average + IM plant in synchronous dq
% Target: 1.5 kW, 4-pole (p=2), 50 Hz base, Vdc-based VSI, closed-loop speed control
% NOTE (important): This step uses actual speed wm for slip computation (debug mode).
% Next step will replace wm with MRAS estimated speed -> sensorless.

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

Tr  = Lr/Rr;        % rotor time constant [s]
Lsigma = Ls - (Lm^2)/Lr;   % leakage (transient) inductance used in current loops

% Inductance matrix (currents -> fluxes)
Lmat = [ Ls   0   Lm   0;
          0  Ls    0  Lm;
         Lm   0   Lr   0;
          0  Lm    0  Lr ];

%% -------------------- DC link voltage and SVPWM limit --------------------
Vdc = 560;                      % [V]
Vmax_ab = Vdc/sqrt(3);          % max |v_alpha_beta| in linear SVPWM region

%% -------------------- Rated load torque baseline -------------------------
P_rated = 1.5e3;
n_sync  = 1500;
w_sync  = 2*pi*n_sync/60;
T_rated = P_rated / w_sync;     % ~9.55 N.m
Tl = T_rated;                   % constant load

%% -------------------- References ----------------------------------------
w_ref_rpm = 1500;               % speed reference [rpm]
w_ref = w_ref_rpm*(2*pi/60);    % [rad/s]

% Flux reference (choose via Id_ref)
% For IM, Id_ref usually near magnetizing current at rated V/f.
% Here choose a reasonable value (tunable):
Id_ref = 5.5;                   % [A] magnetizing current reference (d-axis)

% Current limit (approx). 1.5 kW at 400V has ~3-4 A RMS; peaks can be higher.
Imax = 12;                      % [A] limit for sqrt(id^2+iq^2)

%% -------------------- Controller gains ----------------------------------
% Current loop bandwidth choice
wc_i = 2*pi*400;                % ~400 Hz electrical bandwidth (stable with Ts=100us)

Kp_id = Lsigma*wc_i;
Ki_id = Rs*wc_i;

Kp_iq = Lsigma*wc_i;
Ki_iq = Rs*wc_i;

% Speed loop gains (simple robust tuning)
% We design a PI from speed error -> torque command.
% Then convert torque command -> iq_ref using Te = (3/2)*p*(Lm/Lr)*psi_r*iq
% Approx rotor flux magnitude: psi_r ~= Lm*Id_ref (rough at steady state)
psi_r_ref = Lm*Id_ref;          % [Wb] approximate (good enough for control design)

Kt = (3/2)*p*(Lm/Lr)*psi_r_ref; % torque constant relating iq to torque (approx)

% Speed loop bandwidth ~ 10-20 Hz (rad/s ~ 60-120)
wc_w = 2*pi*12;                 % ~12 Hz speed loop bandwidth
Kp_w = (J*wc_w)/Kt;
Ki_w = (J*wc_w^2)/(Kt*4);       % conservative integral action

% Anti-windup clamps
Tmax = 3*T_rated;               % torque limit [N.m]
Vdq_max = 0.95*Vmax_ab;         % voltage vector magnitude limit in alpha-beta (approx)

%% -------------------- Simulation settings --------------------------------
Ts   = 1e-4;        % 100 us
Tend = 3.0;
t = (0:Ts:Tend).';
N = numel(t);

%% -------------------- States --------------------------------------------
% Plant states: x = [ids iqs idr iqr wm]
x = zeros(5,1);

% Controller integrators
int_w  = 0;
int_id = 0;
int_iq = 0;

% Electrical angle for synchronous reference frame used by controller+plant
theta_e = 0;

%% -------------------- Logs ----------------------------------------------
rpm = zeros(N,1);
rpm_ref = w_ref_rpm*ones(N,1);
Te  = zeros(N,1);

ids = zeros(N,1); iqs = zeros(N,1);
id_ref_log = zeros(N,1); iq_ref_log = zeros(N,1);

vds_log = zeros(N,1); vqs_log = zeros(N,1);
da_log = zeros(N,1); db_log = zeros(N,1); dc_log = zeros(N,1);

%% -------------------- Main loop -----------------------------------------
for k = 1:N
    tk = t(k);

    % Unpack plant states (in synchronous dq frame)
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_k  = x(5);          % mechanical rad/s (ACTUAL, used here in slip calc)

    % ---------------- SPEED CONTROLLER (PI) ----------------
    ew = (w_ref - wm_k);

    % PI -> torque command
    int_w = int_w + Ts*ew;
    Tcmd = Kp_w*ew + Ki_w*int_w;

    % torque clamp + anti-windup (simple back-calculation style)
    if Tcmd > Tmax
        Tcmd = Tmax;
        int_w = (Tcmd - Kp_w*ew)/Ki_w;
    elseif Tcmd < -Tmax
        Tcmd = -Tmax;
        int_w = (Tcmd - Kp_w*ew)/Ki_w;
    end

    % Convert torque command -> iq reference (using approximate torque constant)
    iq_ref = Tcmd / max(Kt, 1e-6);

    % Flux (Id) reference
    id_ref = Id_ref;

    % Current magnitude limiting (keep inside Imax circle)
    Iref_mag = hypot(id_ref, iq_ref);
    if Iref_mag > Imax
        scale = Imax / Iref_mag;
        id_ref = id_ref * scale;
        iq_ref = iq_ref * scale;
    end

    % ---------------- SLIP + ELECTRICAL SPEED ----------------
    % Indirect FOC slip frequency:
    % w_sl = (1/Tr) * (Lm/psi_r_ref) * iq_ref
    % Use psi_r_ref ~ Lm*id_ref (consistent with your chosen magnetizing current)
    psi_r_cmd = max(Lm*abs(id_ref), 0.05);    % avoid divide by zero
    w_sl = (1/Tr)*(Lm/psi_r_cmd)*iq_ref;

    w_r_e = p*wm_k;           % rotor electrical speed
    w_e   = w_r_e + w_sl;     % synchronous electrical speed

    theta_e = wrapTo2Pi(theta_e + Ts*w_e);

    % ---------------- CURRENT CONTROLLERS (PI) ----------------
    eid = id_ref - ids_k;
    eiq = iq_ref - iqs_k;

    int_id = int_id + Ts*eid;
    int_iq = int_iq + Ts*eiq;

    % PI outputs (voltage across Lsigma dynamics)
    v_id_pi = Kp_id*eid + Ki_id*int_id;
    v_iq_pi = Kp_iq*eiq + Ki_iq*int_iq;

    % Decoupling / feedforward (simple, improves stability)
    % Using approximate stator voltage model with leakage inductance:
    % vds ≈ Rs*ids + v_id_pi - w_e*Lsigma*iqs
    % vqs ≈ Rs*iqs + v_iq_pi + w_e*(Lsigma*ids + (Lm/Lr)*psi_r_cmd )
    vds = Rs*ids_k + v_id_pi - w_e*Lsigma*iqs_k;
    vqs = Rs*iqs_k + v_iq_pi + w_e*(Lsigma*ids_k + (Lm/Lr)*psi_r_cmd);

    % Voltage limitation in alpha-beta (SVPWM linear limit)
    [v_alpha, v_beta] = dq_to_alphabeta(vds, vqs, theta_e);
    Vmag = hypot(v_alpha, v_beta);
    if Vmag > Vdq_max
        sc = Vdq_max / Vmag;
        v_alpha = v_alpha*sc;
        v_beta  = v_beta*sc;
        % anti-windup for current integrators: pull back
        vds = vds*sc; vqs = vqs*sc;
        int_id = (vds - Rs*ids_k + w_e*Lsigma*iqs_k - Kp_id*eid)/Ki_id;
        int_iq = (vqs - Rs*iqs_k - w_e*(Lsigma*ids_k + (Lm/Lr)*psi_r_cmd) - Kp_iq*eiq)/Ki_iq;
    end

    % ---------------- SVPWM (average) + voltage reconstruction ------------
    [da, db, dc] = svpwm_duty(v_alpha, v_beta, Vdc);
    [va_k, vb_k, vc_k] = duty_to_vabc(da, db, dc, Vdc);

    % Applied stator voltages in controller synchronous dq frame
    [vds_app, vqs_app] = abc_to_dq(va_k, vb_k, vc_k, theta_e);

    % ---------------- PLANT UPDATE (IM dq in synchronous frame) -----------
    % Fluxes from currents
    i4  = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat*i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);

    % Electromagnetic torque
    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k);

    % Slip for plant rotor equations uses actual rotor speed:
    wsl_plant = w_e - p*wm_k;

    % dq flux derivatives in synchronous frame
    dpsi = zeros(4,1);
    dpsi(1) = vds_app - Rs*ids_k + w_e*psiqs;
    dpsi(2) = vqs_app - Rs*iqs_k - w_e*psids;
    dpsi(3) =         - Rr*idr_k + wsl_plant*psiqr;
    dpsi(4) =         - Rr*iqr_k - wsl_plant*psidr;

    % Convert dpsi -> di
    di = Lmat \ dpsi;

    % Mechanical dynamics
    dwm = (Te_k - Tl - B*wm_k)/J;

    % Integrate
    x(1:4) = x(1:4) + Ts*di;
    x(5)   = x(5)   + Ts*dwm;

    % ---------------- LOGS -----------------------------------------------
    rpm(k) = wm_k*60/(2*pi);
    Te(k)  = Te_k;

    ids(k) = ids_k; iqs(k) = iqs_k;
    id_ref_log(k) = id_ref;
    iq_ref_log(k) = iq_ref;

    vds_log(k) = vds_app;
    vqs_log(k) = vqs_app;

    da_log(k)=da; db_log(k)=db; dc_log(k)=dc;
end

%% -------------------- Plots ---------------------------------------------
figure; plot(t, rpm_ref, 'LineWidth', 1.2); hold on;
plot(t, rpm, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title('FOC Speed Control (Healthy, SVPWM avg)');
legend('Reference','Actual','Location','best');

figure; plot(t, Te, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Torque [N·m]');
title('Electromagnetic Torque');

figure; plot(t, ids, 'LineWidth', 1.2); hold on;
plot(t, iqs, 'LineWidth', 1.2); grid on;
plot(t, id_ref_log, '--', 'LineWidth', 1.0);
plot(t, iq_ref_log, '--', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Current [A]');
title('Stator dq currents (FOC)');
legend('i_{ds}','i_{qs}','i_{ds}^*','i_{qs}^*','Location','best');

figure;
idx = t<=0.04;
plot(t(idx), da_log(idx), 'LineWidth', 1.0); hold on;
plot(t(idx), db_log(idx), 'LineWidth', 1.0);
plot(t(idx), dc_log(idx), 'LineWidth', 1.0);
grid on; xlabel('Time [s]'); ylabel('Duty');
title('SVPWM duty ratios (first 40 ms)');
legend('d_a','d_b','d_c','Location','best');

%% -------------------- Local functions -----------------------------------
function [da, db, dc] = svpwm_duty(v_alpha, v_beta, Vdc)
% SVPWM duty computation (average model, 2-level VSI).
    Vmax = Vdc/sqrt(3);
    Vref = hypot(v_alpha, v_beta);
    if Vref > Vmax
        s = Vmax / Vref;
        v_alpha = v_alpha*s;
        v_beta  = v_beta*s;
    end

    % Convert alpha-beta voltage to equivalent normalized 3-phase commands
    ua = (2/3)*( v_alpha )/Vdc;
    ub = (2/3)*(-v_alpha/2 + (sqrt(3)/2)*v_beta)/Vdc;
    uc = (2/3)*(-v_alpha/2 - (sqrt(3)/2)*v_beta)/Vdc;

    u_max = max([ua, ub, uc]);
    u_min = min([ua, ub, uc]);
    u0 = -0.5*(u_max + u_min);    % zero-sequence injection

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

function [v_alpha, v_beta] = dq_to_alphabeta(vd, vq, theta)
% inverse Park: dq -> alpha-beta
    c = cos(theta); s = sin(theta);
    v_alpha =  vd*c - vq*s;
    v_beta  =  vd*s + vq*c;
end
