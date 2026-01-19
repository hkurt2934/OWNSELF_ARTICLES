%% IM_Sensorless_FOC_FTC_WholeLeg.m
% Sensorless Induction Motor Speed Control:
%   - Healthy: Indirect Rotor-Flux-Oriented Control (IRFOC) + MRAS speed estimator
%   - Fault (whole-leg open, i.e., two-switch fault on phase-A leg): FTC mode
%       => v_alpha = 0, ia = 0, only v_beta (via v_bc) is effective + flux constraint
%
% Plant model is the same style as your file (complex alpha-beta flux states).

clear; clc; close all;

%% ================== Machine Parameters ==================
Rs = 0.435;
Rr = 0.816;
Lls = 0.002;
Llr = 0.002;
Lm  = 0.0693;
Ls  = Lls + Lm;
Lr  = Llr + Lm;
P   = 4;
J   = 0.089;
B   = 0.001;

Delta = Ls*Lr - Lm^2;
Lsigma = Ls - (Lm^2)/Lr;        % leakage inductance (sigma)
Tr = Lr/Rr;                     % rotor time constant

%% ================== Simulation Settings ==================
Ts    = 1e-4;
Tend  = 3.0;
time  = 0:Ts:Tend;
N     = numel(time);

tfault = 1.0;                   % fault time
TL     = 50;                    % load torque [N.m]

w_ref_rpm = 1500;               % speed reference [rpm]
w_ref     = w_ref_rpm*2*pi/60;  % [rad/s] mechanical

%% ================== Inverter Limit ==================
Vdc     = 600;
Vpk_max = Vdc/sqrt(3);          % fundamental phase peak limit (SVPWM approx)

%% ================== Control Targets ==================
psi_r_ref = 0.90;               % rotor flux reference [Wb] (tune)
Id_ref    = psi_r_ref / Lm;     % i_d* ~ psi_r/Lm

Iq_max    = 200;                % torque current limit [A]
Id_max    = 50;                 % flux current limit  [A]

%% ================== Current PI (dq) Tuning ==================
% Simple bandwidth-based choice:
%   Kp = Lsigma*wci, Ki = Rs*wci
fci = 500;                      % current loop bandwidth [Hz] (tune)
wci = 2*pi*fci;

Kp_i = Lsigma*wci;
Ki_i = Rs*wci;

%% ================== Speed PI Tuning ==================
% i_q_ref from speed error (sensorless speed estimate)
% Practical, saturated PI (tune as needed)
Kp_w = 3.4;      % for ~8 Hz-like response with your parameters (tune)
Ki_w = 85.7;     % (tune)

%% ================== MRAS Speed Estimator ==================
% Voltage model (stator flux) with "leaky integrator" to reduce drift
Tf_vm = 0.5;                     % leak time constant [s] (tune)

% MRAS adaptation (PI)
Kp_mras = 30;                    % (tune)
Ki_mras = 2000;                  % (tune)
omega_r_elec_max = 2*pi*300;     % electrical rad/s clamp

%% ================== Run 3 cases ==================
% 1) Healthy (no fault)
% 2) Fault without FTC (controller not reconfigured)
% 3) Fault with FTC (reconfigured single-axis control)
res_h  = sim_case(false,false);
res_fw = sim_case(true ,false);
res_ft = sim_case(true ,true );

%% ================== Quick Metrics ==================
win = time > 2.0;

fprintf('\n=== Healthy (t>2s) ===\n');
print_metrics(res_h, win);

fprintf('\n=== Fault WITHOUT FTC (t>2s) ===\n');
print_metrics(res_fw, win);

fprintf('\n=== Fault WITH FTC (t>2s) ===\n');
print_metrics(res_ft, win);

%% ================== Plots ==================
figure('Name','Speed');
plot(time,res_h.wm_rpm,'LineWidth',1.1); hold on;
plot(time,res_fw.wm_rpm,'LineWidth',1.1);
plot(time,res_ft.wm_rpm,'LineWidth',1.1);
xline(tfault,'r--','Fault'); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title('Rotor Speed: Healthy vs Fault (w/o FTC) vs Fault (with FTC)');
legend('Healthy','Fault w/o FTC','Fault with FTC','Location','best');

figure('Name','Torque');
plot(time,res_h.Te,'LineWidth',1.0); hold on;
plot(time,res_fw.Te,'LineWidth',1.0);
plot(time,res_ft.Te,'LineWidth',1.0);
yline(TL,'k--','Load'); xline(tfault,'r--','Fault'); grid on;
xlabel('Time [s]'); ylabel('Torque [N.m]');
title('Electromagnetic Torque');
legend('Healthy','Fault w/o FTC','Fault with FTC','Location','best');

figure('Name','Phase Currents (Fault with FTC)');
subplot(3,1,1); plot(time,res_ft.ia,'LineWidth',0.9); xline(tfault,'r--'); grid on; ylabel('i_a [A]'); title('Fault with FTC: i_a,i_b,i_c');
subplot(3,1,2); plot(time,res_ft.ib,'LineWidth',0.9); xline(tfault,'r--'); grid on; ylabel('i_b [A]');
subplot(3,1,3); plot(time,res_ft.ic,'LineWidth',0.9); xline(tfault,'r--'); grid on; ylabel('i_c [A]'); xlabel('Time [s]');

figure('Name','Sensorless Speed Estimation (Fault with FTC)');
plot(time,res_ft.wm_rpm,'LineWidth',1.1); hold on;
plot(time,res_ft.wm_hat_rpm,'LineWidth',1.1);
xline(tfault,'r--','Fault'); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title('Measured (plant) vs Estimated Speed (MRAS)');
legend('Plant speed','Estimated speed','Location','best');

%% ================== Local Functions ==================

function out = sim_case(enable_fault, enable_FTC)

    % --- pull base variables
    Rs = evalin('base','Rs'); Rr = evalin('base','Rr');
    Lm = evalin('base','Lm'); Ls = evalin('base','Ls'); Lr = evalin('base','Lr');
    Delta = evalin('base','Delta'); Lsigma = evalin('base','Lsigma'); Tr = evalin('base','Tr');
    P  = evalin('base','P');  J  = evalin('base','J'); B = evalin('base','B');

    Ts = evalin('base','Ts'); time = evalin('base','time'); N = evalin('base','N');
    tfault = evalin('base','tfault'); TL = evalin('base','TL'); w_ref = evalin('base','w_ref');

    Vpk_max = evalin('base','Vpk_max');

    psi_r_ref = evalin('base','psi_r_ref');
    Id_ref    = evalin('base','Id_ref');
    Iq_max    = evalin('base','Iq_max');
    Id_max    = evalin('base','Id_max');

    Kp_i = evalin('base','Kp_i'); Ki_i = evalin('base','Ki_i');
    Kp_w = evalin('base','Kp_w'); Ki_w = evalin('base','Ki_w');

    Tf_vm = evalin('base','Tf_vm');
    Kp_mras = evalin('base','Kp_mras'); Ki_mras = evalin('base','Ki_mras');
    omega_r_elec_max = evalin('base','omega_r_elec_max');

    % --- plant states (complex alpha-beta)
    % initialize with some flux to help sensorless start
    psi_r = psi_r_ref + 1j*0;
    psi_s = (Lm/Lr)*real(psi_r) + 1j*0;
    wm    = 0;

    % --- controller integrators
    int_w   = 0;
    int_id  = 0; int_iq  = 0;
    int_ib  = 0; % for FTC single-axis beta current PI

    % --- estimator states
    psi_s_vm  = 0 + 1j*0;     % stator flux (voltage model)
    psi_r_cm  = psi_r;        % rotor flux (current model)
    omega_r_hat = 0;          % estimated rotor electrical speed [rad/s]
    int_mras    = 0;
    theta_e     = 0;          % synchronous electrical angle

    % --- logs
    ia = zeros(1,N); ib = zeros(1,N); ic = zeros(1,N);
    wm_rpm = zeros(1,N); wm_hat_rpm = zeros(1,N);
    Te_log = zeros(1,N);

    for k = 1:N
        t = time(k);

        fault = enable_fault && (t >= tfault);
        use_FTC = fault && enable_FTC;

        omega_r = (P/2)*wm;  % plant electrical speed

        % ---------- fault constraint (plant) before step ----------
        if fault
            % same constraint idea as your file (enforces i_alpha ~ 0)
            psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
        end

        % ---------- measure currents from plant states ----------
        i_s = (Lr*psi_s - Lm*psi_r)/Delta;
        i_alpha = real(i_s);
        i_beta  = imag(i_s);

        % phase currents (for plotting)
        if fault
            ia(k) = 0;
            ib(k) =  (sqrt(3)/2)*i_beta;
            ic(k) = -(sqrt(3)/2)*i_beta;
        else
            [ia(k), ib(k), ic(k)] = inv_clarke(i_alpha, i_beta);
        end

        % ---------- Sensorless MRAS ----------
        % We'll use the *applied* v_s (ideal inverter assumption) for VM.
        % Compute control first, then update estimator.

        % ---------- Outer speed loop (uses estimated speed) ----------
        wm_hat = omega_r_hat/(P/2);     % mechanical rad/s
        e_w = w_ref - wm_hat;
        int_w = int_w + e_w*Ts;

        i_q_ref = Kp_w*e_w + Ki_w*int_w;
        i_q_ref = sat(i_q_ref, Iq_max);

        % anti-windup (simple)
        if abs(i_q_ref) >= Iq_max
            int_w = int_w - ( (Kp_w*e_w + Ki_w*int_w) - i_q_ref )/max(Ki_w,1e-9);
        end

        i_d_ref = sat(Id_ref, Id_max);

        % ---------- Park transform using theta_e ----------
        [i_d, i_q] = park(i_alpha, i_beta, theta_e);

        % ---------- Slip + synchronous speed estimate (IRFOC) ----------
        psi_r_mag = max(abs(psi_r_cm), 1e-6);
        % rotor-flux-oriented: omega_slip = (Lm/Tr)*(i_q/psi_r_d)
        % Here psi_r_d ~ |psi_r|
        omega_slip = (Lm/Tr) * (i_q / psi_r_mag);
        omega_syn  = omega_r_hat + omega_slip;

        % ---------- Current control ----------
        if ~use_FTC
            % dq current PI + simple decoupling/feedforward
            e_id = i_d_ref - i_d;
            e_iq = i_q_ref - i_q;

            int_id = int_id + e_id*Ki_i*Ts;
            int_iq = int_iq + e_iq*Ki_i*Ts;

            v_d_pi = Kp_i*e_id + int_id;
            v_q_pi = Kp_i*e_iq + int_iq;

            % decoupling (approx)
            v_d = v_d_pi - omega_syn*Lsigma*i_q;
            v_q = v_q_pi + omega_syn*(Lsigma*i_d + (Lm/Lr)*psi_r_mag);

            % voltage limit (alpha-beta magnitude)
            [v_alpha, v_beta] = inv_park(v_d, v_q, theta_e);
            [v_alpha, v_beta] = limit_v(v_alpha, v_beta, Vpk_max);

            v_s = v_alpha + 1j*v_beta;

        else
            % ---------- FTC mode (whole-leg open): single-axis beta control ----------
            % We can only command v_beta effectively (v_alpha=0).
            % Use i_beta as the controlled variable; map speed PI torque demand to i_beta_ref.
            i_beta_ref = sat(i_q_ref, Iq_max); % treat i_beta ~ torque current

            e_ib = i_beta_ref - i_beta;
            int_ib = int_ib + e_ib*Ki_i*Ts;
            v_beta_cmd = Kp_i*e_ib + int_ib;

            % limit v_beta to inverter capability
            v_beta_cmd = sat(v_beta_cmd, Vpk_max);

            v_s = 1j*v_beta_cmd; % v_alpha = 0
        end

        % ---------- Voltage model stator flux (leaky) ----------
        % psi_s_vm_dot = v_s - Rs*i_s - psi_s_vm/Tf_vm
        psi_s_vm = psi_s_vm + ( (v_s - Rs*i_s) - psi_s_vm/max(Tf_vm,1e-6) )*Ts;

        % rotor flux from voltage model:
        % psi_r_vm = (Lr/Lm)*(psi_s_vm - Lsigma*i_s)
        psi_r_vm = (Lr/Lm)*(psi_s_vm - Lsigma*i_s);

        % current model rotor flux update (uses omega_r_hat)
        dpsi_r_cm = (Lm/Tr)*i_s - (1/Tr)*psi_r_cm + 1j*omega_r_hat*psi_r_cm;
        psi_r_cm = psi_r_cm + dpsi_r_cm*Ts;

        % MRAS error: cross product between psi_r_vm and psi_r_cm
        e_mras = imag(conj(psi_r_vm)*psi_r_cm);
        int_mras = int_mras + e_mras*Ts;

        omega_r_hat = omega_r_hat + (Kp_mras*e_mras + Ki_mras*int_mras)*Ts;
        omega_r_hat = sat(omega_r_hat, omega_r_elec_max);

        % update synchronous angle
        theta_e = wrap_pi(theta_e + omega_syn*Ts);

        % ---------- Plant flux dynamics ----------
        dpsi_s = v_s - Rs*i_s;
        i_r = (-Lm*psi_s + Ls*psi_r)/Delta;
        dpsi_r = -Rr*i_r + 1j*omega_r*psi_r;

        psi_s = psi_s + dpsi_s*Ts;
        psi_r = psi_r + dpsi_r*Ts;

        % re-apply constraint after integration
        if fault
            psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
        end

        % ---------- torque + mechanical ----------
        i_s = (Lr*psi_s - Lm*psi_r)/Delta;
        Te = 1.5*(P/2) * imag(conj(psi_s)*i_s);
        wm = wm + (Te - TL - B*wm)/J * Ts;

        wm_rpm(k)     = wm*60/(2*pi);
        wm_hat_rpm(k) = (omega_r_hat/(P/2))*60/(2*pi);
        Te_log(k) = Te;
    end

    out.wm_rpm = wm_rpm;
    out.wm_hat_rpm = wm_hat_rpm;
    out.Te = Te_log;
    out.ia = ia; out.ib = ib; out.ic = ic;
end

function print_metrics(res, win)
    Te_mean = mean(res.Te(win));
    Te_pp   = max(res.Te(win)) - min(res.Te(win));
    w_mean  = mean(res.wm_rpm(win));
    w_pp    = max(res.wm_rpm(win)) - min(res.wm_rpm(win));
    ia_rms  = rms(res.ia(win));
    ib_rms  = rms(res.ib(win));
    ic_rms  = rms(res.ic(win));

    fprintf("Te_mean=%.2f N.m, Te_pp=%.2f N.m\n", Te_mean, Te_pp);
    fprintf("w_mean=%.2f rpm, w_ripple_pp=%.2f rpm\n", w_mean, w_pp);
    fprintf("RMS currents: ia=%.2f A, ib=%.2f A, ic=%.2f A\n", ia_rms, ib_rms, ic_rms);
end

%% ====== Utility: Clarke / inverse Clarke ======
function [i_alpha, i_beta] = clarke(ia, ib, ic)
    % power-invariant Clarke
    i_alpha = (2/3)*(ia - 0.5*ib - 0.5*ic);
    i_beta  = (2/3)*((sqrt(3)/2)*(ib - ic));
end

function [ia, ib, ic] = inv_clarke(i_alpha, i_beta)
    ia = i_alpha;
    ib = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
    ic = -0.5*i_alpha - (sqrt(3)/2)*i_beta;
end

%% ====== Park / inverse Park ======
function [d,q] = park(alpha, beta, theta)
    c = cos(theta); s = sin(theta);
    d =  c*alpha + s*beta;
    q = -s*alpha + c*beta;
end

function [alpha,beta] = inv_park(d,q,theta)
    c = cos(theta); s = sin(theta);
    alpha =  c*d - s*q;
    beta  =  s*d + c*q;
end

%% ====== Voltage magnitude limiter in alpha-beta ======
function [v_alpha, v_beta] = limit_v(v_alpha, v_beta, Vmax)
    mag = hypot(v_alpha, v_beta);
    if mag > Vmax
        sc = Vmax/mag;
        v_alpha = v_alpha*sc;
        v_beta  = v_beta*sc;
    end
end

%% ====== Saturation and angle wrap ======
function y = sat(x, lim)
    y = min(max(x, -lim), lim);
end

function th = wrap_pi(th)
    th = mod(th + pi, 2*pi) - pi;
end
