%% STEP 5 - Healthy IM: FOC (IFOC) + SVPWM(avg) + MRAS speed estimator (sensorless)
% Only measured signals assumed: ia, ib, ic and Vdc (and switching duty).
% In simulation: vabc is reconstructed from duties+Vdc (as in real drive).
% Speed loop feedback uses wm_hat from MRAS (NO actual speed used in control).

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
Tr  = Lr/Rr;                       % rotor time constant
Lsigma = Ls - (Lm^2)/Lr;           % transient/leakage inductance

Lmat = [ Ls   0   Lm   0;
          0  Ls    0  Lm;
         Lm   0   Lr   0;
          0  Lm    0  Lr ];

%% -------------------- DC link + SVPWM limits ----------------------------
Vdc = 560;                      % [V]
Vmax_ab = Vdc/sqrt(3);          % SVPWM linear limit for |v_alpha_beta|
Vdq_max = 0.95*Vmax_ab;

%% -------------------- Load torque ---------------------------------------
P_rated = 1.5e3;
n_sync  = 1500;
w_sync  = 2*pi*n_sync/60;
T_rated = P_rated / w_sync;     % ~9.55 N.m
Tl = T_rated;

%% -------------------- References ----------------------------------------
w_ref_rpm = 1500;
w_ref = w_ref_rpm*(2*pi/60);

Id_ref = 4.5;                   % [A] (flux current ref) tuneable
Imax   = 12;                    % [A] current limit
Tmax   = 3*T_rated;             % torque limit

%% -------------------- Current loop PI gains -----------------------------
wc_i = 2*pi*400;
Kp_id = Lsigma*wc_i;  Ki_id = Rs*wc_i;
Kp_iq = Lsigma*wc_i;  Ki_iq = Rs*wc_i;

%% -------------------- Speed loop PI gains (uses wm_hat!) ----------------
% Conservative tuning
psi_r_ref = max(Lm*Id_ref, 0.05);
Kt = (3/2)*p*(Lm/Lr)*psi_r_ref;
wc_w = 2*pi*12;
Kp_w = (J*wc_w)/Kt;
Ki_w = (J*wc_w^2)/(Kt*4);

%% -------------------- MRAS gains (tunable) -------------------------------
% These work as a starting point for Ts=100us and this motor set.
Kp_mras = 80;           % speed adaptation proportional gain
Ki_mras = 3000;         % speed adaptation integral gain
w_r_hat_lim = 2*pi*3000/60 * p;   % electrical rad/s limit (~3000 rpm mech)

% Voltage-model flux anti-drift leakage (helps long simulations)
Tf_psi = 1.0;           % [s] leakage time constant
alpha_psi = 1/Tf_psi;

%% -------------------- Simulation settings --------------------------------
Ts   = 1e-4;
Tend = 3.0;
t = (0:Ts:Tend).';
N = numel(t);

%% -------------------- Plant states (dq synchronous) ----------------------
% x = [ids iqs idr iqr wm_actual]
x = zeros(5,1);

%% -------------------- Controller integrators -----------------------------
int_w  = 0;
int_id = 0;
int_iq = 0;

%% -------------------- Electrical angle used by control -------------------
theta_e = 0;

%% -------------------- MRAS states (stationary alpha-beta) ----------------
psi_s_ab   = [0;0];    % stator flux (voltage model)
psi_rv_ab  = [0;0];    % rotor flux from voltage model
psi_ri_ab  = [0;0];    % rotor flux from current model
int_mras   = 0;        % integrator for speed adaptation
w_r_hat    = 0;        % estimated rotor electrical speed (rad/s)
wm_hat     = 0;        % estimated mechanical speed (rad/s)

%% -------------------- Logs ----------------------------------------------
rpm_act = zeros(N,1);
rpm_hat = zeros(N,1);
rpm_ref = w_ref_rpm*ones(N,1);

Te  = zeros(N,1);
ids = zeros(N,1); iqs = zeros(N,1);
idref_log = zeros(N,1); iqref_log = zeros(N,1);
da_log = zeros(N,1); db_log = zeros(N,1); dc_log = zeros(N,1);
mras_err = zeros(N,1);

%% -------------------- Main loop -----------------------------------------
for k = 1:N
    tk = t(k);

    % -------- Unpack plant states (dq synchronous frame) ----------
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_act_k  = x(5);  % only for plant mech update, NOT used by control

    % -------- Speed PI uses estimated speed (sensorless) ----------
    ew = (w_ref - wm_hat);

    % PI torque command (with saturation + conditional anti-windup)
    Tcmd_unsat = Kp_w*ew + Ki_w*int_w;
    Tcmd = min(max(Tcmd_unsat, -Tmax), Tmax);

    isSatHigh = (Tcmd_unsat >  Tmax);
    isSatLow  = (Tcmd_unsat < -Tmax);
    if (~isSatHigh && ~isSatLow)
        int_w = int_w + Ts*ew;
    else
        if (isSatHigh && ew < 0), int_w = int_w + Ts*ew; end
        if (isSatLow  && ew > 0), int_w = int_w + Ts*ew; end
    end

    % torque -> iq_ref (use current rotor flux command)
    psi_r_cmd = max(Lm*abs(Id_ref), 0.05);
    Kt_inst = (3/2)*p*(Lm/Lr)*psi_r_cmd;
    iq_ref = Tcmd / max(Kt_inst, 1e-6);
    id_ref = Id_ref;

    % current magnitude limiting
    Iref_mag = hypot(id_ref, iq_ref);
    if Iref_mag > Imax
        sc = Imax/Iref_mag;
        id_ref = id_ref*sc;
        iq_ref = iq_ref*sc;
    end

    % -------- Indirect FOC slip using estimated speed/flux --------
    % slip frequency based on commanded flux (Id_ref) and iq_ref
    psi_r_for_slip = max(Lm*abs(id_ref), 0.05);
    w_sl = (1/Tr)*(Lm/psi_r_for_slip)*iq_ref;

    % synchronous electrical speed uses estimated rotor electrical speed
    w_e = w_r_hat + w_sl;
    theta_e = wrapTo2Pi(theta_e + Ts*w_e);

    % -------- Current controllers (dq) ----------
    eid = id_ref - ids_k;
    eiq = iq_ref - iqs_k;
    int_id = int_id + Ts*eid;
    int_iq = int_iq + Ts*eiq;

    v_id_pi = Kp_id*eid + Ki_id*int_id;
    v_iq_pi = Kp_iq*eiq + Ki_iq*int_iq;

    % decoupling
    vds = Rs*ids_k + v_id_pi - w_e*Lsigma*iqs_k;
    vqs = Rs*iqs_k + v_iq_pi + w_e*(Lsigma*ids_k + (Lm/Lr)*psi_r_for_slip);

    % dq -> alpha-beta for SVPWM limiting
    [v_alpha_cmd, v_beta_cmd] = dq_to_alphabeta(vds, vqs, theta_e);
    Vmag = hypot(v_alpha_cmd, v_beta_cmd);
    if Vmag > Vdq_max
        sc = Vdq_max / Vmag;
        v_alpha_cmd = v_alpha_cmd*sc;
        v_beta_cmd  = v_beta_cmd*sc;
    end

    % -------- SVPWM duties + reconstruct vabc from (duty,Vdc) ----------
    [da, db, dc] = svpwm_duty(v_alpha_cmd, v_beta_cmd, Vdc);
    [va, vb, vc] = duty_to_vabc(da, db, dc, Vdc);

    % apply to plant: abc -> dq with theta_e (control frame)
    [vds_app, vqs_app] = abc_to_dq(va, vb, vc, theta_e);

    % -------- Plant update: IM model in synchronous dq ----------
    i4 = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat*i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);

    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k);

    wsl_plant = w_e - p*wm_act_k;

    dpsi = zeros(4,1);
    dpsi(1) = vds_app - Rs*ids_k + w_e*psiqs;
    dpsi(2) = vqs_app - Rs*iqs_k - w_e*psids;
    dpsi(3) =         - Rr*idr_k + wsl_plant*psiqr;
    dpsi(4) =         - Rr*iqr_k - wsl_plant*psidr;

    di = Lmat \ dpsi;
    dwm = (Te_k - Tl - B*wm_act_k)/J;

    x(1:4) = x(1:4) + Ts*di;
    x(5)   = x(5)   + Ts*dwm;

    % ============================================================
    % ================= MRAS SPEED ESTIMATOR ======================
    % ============================================================
    % Measurements used: vabc reconstructed (va,vb,vc), and iabc measured.
    % In simulation iabc is reconstructed from ids/iqs + theta_e.

    % Reconstruct iabc (measurement) from dq currents
    [ia, ib, ic] = dq_to_abc(ids_k, iqs_k, theta_e);

    % Clarke to alpha-beta
    [v_alpha, v_beta] = abc_to_alphabeta(va, vb, vc);
    [i_alpha, i_beta] = abc_to_alphabeta(ia, ib, ic);
    i_ab = [i_alpha; i_beta];
    v_ab = [v_alpha; v_beta];

    % Voltage-model stator flux (with small leakage to avoid drift)
    % dpsi_s = v - Rs*i - alpha_psi*psi_s
    dpsi_s = v_ab - Rs*i_ab - alpha_psi*psi_s_ab;
    psi_s_ab = psi_s_ab + Ts*dpsi_s;

    % Rotor flux from voltage model
    psi_rv_ab = (Lr/Lm)*(psi_s_ab - Lsigma*i_ab);

    % Current-model rotor flux (adjustable model)
    % dpsi_ri = (Lm/Tr)*i - (1/Tr)*psi_ri + w_r_hat * J * psi_ri
    Jmat = [0 -1; 1 0];
    dpsi_ri = (Lm/Tr)*i_ab - (1/Tr)*psi_ri_ab + w_r_hat*(Jmat*psi_ri_ab);
    psi_ri_ab = psi_ri_ab + Ts*dpsi_ri;

    % MRAS error (cross product)
    e_mras = psi_rv_ab(1)*psi_ri_ab(2) - psi_rv_ab(2)*psi_ri_ab(1);

    % PI adaptation for rotor electrical speed estimate
    int_mras = int_mras + Ts*e_mras;
    w_r_hat = Kp_mras*e_mras + Ki_mras*int_mras;

    % clamp estimated speed
    w_r_hat = min(max(w_r_hat, -w_r_hat_lim), w_r_hat_lim);

    % mechanical speed estimate
    wm_hat = w_r_hat / p;

    % -------- Logs ----------
    rpm_act(k) = wm_act_k*60/(2*pi);
    rpm_hat(k) = wm_hat*60/(2*pi);
    Te(k) = Te_k;
    ids(k)=ids_k; iqs(k)=iqs_k;
    idref_log(k)=id_ref; iqref_log(k)=iq_ref;
    da_log(k)=da; db_log(k)=db; dc_log(k)=dc;
    mras_err(k)=e_mras;
end

%% -------------------- Plots ---------------------------------------------
figure; plot(t, rpm_ref, 'LineWidth', 1.2); hold on;
plot(t, rpm_act, 'LineWidth', 1.2);
plot(t, rpm_hat, '--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Speed [rpm]');
title('Sensorless FOC Speed (MRAS) - Reference vs Actual vs Estimated');
legend('Reference','Actual','Estimated','Location','best');

figure; plot(t, Te, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Torque [N·m]');
title('Electromagnetic Torque');

figure; plot(t, ids, 'LineWidth', 1.2); hold on;
plot(t, iqs, 'LineWidth', 1.2);
plot(t, idref_log, '--', 'LineWidth', 1.0);
plot(t, iqref_log, '--', 'LineWidth', 1.0);
grid on; xlabel('Time [s]'); ylabel('Current [A]');
title('Stator dq currents (Sensorless FOC)');
legend('i_{ds}','i_{qs}','i_{ds}^*','i_{qs}^*','Location','best');

figure; plot(t, mras_err, 'LineWidth', 1.1); grid on;
xlabel('Time [s]'); ylabel('MRAS error (cross-product)');
title('MRAS adaptation error');

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
    Vmax = Vdc/sqrt(3);
    Vref = hypot(v_alpha, v_beta);
    if Vref > Vmax
        s = Vmax / Vref;
        v_alpha = v_alpha*s;
        v_beta  = v_beta*s;
    end

    ua = (2/3)*( v_alpha )/Vdc;
    ub = (2/3)*(-v_alpha/2 + (sqrt(3)/2)*v_beta)/Vdc;
    uc = (2/3)*(-v_alpha/2 - (sqrt(3)/2)*v_beta)/Vdc;

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
    v_aO = (2*da - 1)*(Vdc/2);
    v_bO = (2*db - 1)*(Vdc/2);
    v_cO = (2*dc - 1)*(Vdc/2);
    v_cm = (v_aO + v_bO + v_cO)/3;
    va = v_aO - v_cm;
    vb = v_bO - v_cm;
    vc = v_cO - v_cm;
end

function [vd, vq] = abc_to_dq(va, vb, vc, theta)
    alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
    c = cos(theta); s = sin(theta);
    vd =  alpha*c + beta*s;
    vq = -alpha*s + beta*c;
end

function [v_alpha, v_beta] = dq_to_alphabeta(vd, vq, theta)
    c = cos(theta); s = sin(theta);
    v_alpha =  vd*c - vq*s;
    v_beta  =  vd*s + vq*c;
end

function [alpha, beta] = abc_to_alphabeta(a, b, c)
    alpha = (2/3)*(a - 0.5*b - 0.5*c);
    beta  = (2/3)*((sqrt(3)/2)*(b - c));
end

function [a, b, c] = dq_to_abc(id, iq, theta)
    % dq -> alpha-beta
    cth = cos(theta); sth = sin(theta);
    i_alpha =  id*cth - iq*sth;
    i_beta  =  id*sth + iq*cth;

    % alpha-beta -> abc (inverse Clarke)
    a = i_alpha;
    b = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
    c = -0.5*i_alpha - (sqrt(3)/2)*i_beta;
end
