%% Test_X1_sensorless_FTC_MRAS.m
% Sensorless Fault-Tolerant Control (FTC) demo for Induction Motor
% using a simple MRAS speed estimator + fault detection + post-fault
% reconfiguration for a whole-leg open-circuit (phase-A leg open).
%
% Scenarios (same load and speed reference):
%   1) Healthy (sensorless)
%   2) Whole-leg OCF WITHOUT FTC (sensorless controller not reconfigured)
%   3) Whole-leg OCF WITH FTC (FDI + reconfiguration)
%
% User case:
%   w_ref = 1500 rpm
%   TL    = 50 N.m
%   tfault = 1 s
%
% Notes:
% - Plant is an averaged flux-linkage IM model (no PWM switching ripple).
% - Whole-leg fault is enforced with ia=0 constraint (ib=-ic), consistent
%   with prior scripts.
%
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
Tr = Lr/Rr;                 % rotor time constant

%% ================== Simulation Settings ==================
Ts = 1e-4;
Tend = 3.0;
time = 0:Ts:Tend;
N = numel(time);

tfault = 1.0;
TL = 50;                      % N.m
w_ref = 1500 * 2*pi/60;       % mechanical rad/s
ws_ref = (P/2) * w_ref;       % electrical rad/s

%% ================== Inverter + Scalar Control ==================
Vdc = 600;
Vpk_max = Vdc/sqrt(3);        % SVPWM fundamental phase peak limit

% V/f base
f_base   = 50;
ws_base  = 2*pi*f_base;
Vpk_base = 325;               % ~230 Vrms phase => 325 Vpk
Kvf      = Vpk_base/ws_base;

% Speed PI (uses estimated speed)
Kpw = 0.6;
Kiw = 30.0;
ws_min = 0;
ws_max = 2*pi*80;

% Post-fault reconfiguration gains (FTC mode)
Kvf_fault = 1.10*Kvf;         % slightly boost V/f after fault
Kiw_fault = 18.0;             % reduce integral action to avoid windup
Kpw_fault = 0.6;              % keep proportional

%% ================== Sensorless MRAS Parameters ==================
% Voltage-model flux integrator leakage to prevent drift
kf_flux = 5.0;                % rad/s (small)

% MRAS adaptation gains (electrical rotor speed)
% Error is normalized => gains can be moderate.
Kp_mras = 150.0;
Ki_mras = 8000.0;
wr_hat_min = 0;
wr_hat_max = ws_max;

%% ================== Fault Detection (FDI) ==================
Twin = 0.02;                         % 20 ms window
Nw = max(20, round(Twin/Ts));
thr_ia  = 0.10;                      % ia_rms/ib_rms threshold
thr_sum = 0.10;                      % RMS(ib+ic)/ib_rms threshold
hold_windows = 3;                    % require N consecutive windows

%% ================== Run scenarios ==================
opt1 = make_opts('Healthy (sensorless)', false, false);
opt2 = make_opts('Fault w/o FTC (sensorless)', true,  false);
opt3 = make_opts('Fault WITH FTC (sensorless)', true,  true);

res1 = sim_case(opt1);
res2 = sim_case(opt2);
res3 = sim_case(opt3);

%% ================== Plots ==================
% Speed (true vs estimated)
figure('Name','Speed - true vs estimated');
subplot(2,1,1);
plot(time, res1.wm_rpm,'LineWidth',1.1); hold on;
plot(time, res2.wm_rpm,'LineWidth',1.1);
plot(time, res3.wm_rpm,'LineWidth',1.1);
xline(tfault,'r--','Fault'); grid on;
title('Rotor Speed (true)'); xlabel('Time [s]'); ylabel('Speed [rpm]');
legend('Healthy','Fault w/o FTC','Fault WITH FTC','Location','best');

subplot(2,1,2);
plot(time, res1.wm_hat_rpm,'LineWidth',1.1); hold on;
plot(time, res2.wm_hat_rpm,'LineWidth',1.1);
plot(time, res3.wm_hat_rpm,'LineWidth',1.1);
xline(tfault,'r--','Fault'); grid on;
title('Estimated Speed (MRAS)'); xlabel('Time [s]'); ylabel('Speed est [rpm]');
legend('Healthy','Fault w/o FTC','Fault WITH FTC','Location','best');

% Torque
figure('Name','Torque');
plot(time, res1.Te,'LineWidth',1.0); hold on;
plot(time, res2.Te,'LineWidth',1.0);
plot(time, res3.Te,'LineWidth',1.0);
xline(tfault,'r--','Fault'); yline(TL,'k--','Load'); grid on;
title('Electromagnetic Torque'); xlabel('Time [s]'); ylabel('Torque [N.m]');
legend('Healthy','Fault w/o FTC','Fault WITH FTC','Location','best');

% Phase currents (overlay)
figure('Name','Phase currents');
subplot(3,1,1);
plot(time, res1.ia,'LineWidth',0.9); hold on;
plot(time, res2.ia,'LineWidth',0.9);
plot(time, res3.ia,'LineWidth',0.9);
xline(tfault,'r--','Fault'); grid on; ylabel('i_a [A]');
legend('Healthy','Fault w/o FTC','Fault WITH FTC','Location','best');
title('Stator Phase Currents');

subplot(3,1,2);
plot(time, res1.ib,'LineWidth',0.9); hold on;
plot(time, res2.ib,'LineWidth',0.9);
plot(time, res3.ib,'LineWidth',0.9);
xline(tfault,'r--','Fault'); grid on; ylabel('i_b [A]');

subplot(3,1,3);
plot(time, res1.ic,'LineWidth',0.9); hold on;
plot(time, res2.ic,'LineWidth',0.9);
plot(time, res3.ic,'LineWidth',0.9);
xline(tfault,'r--','Fault'); grid on; ylabel('i_c [A]'); xlabel('Time [s]');

% Fault flag
figure('Name','FTC supervisor');
stairs(time, res3.ftc_active,'LineWidth',1.2); hold on;
xline(tfault,'r--','Fault'); grid on;
ylabel('FTC active'); xlabel('Time [s]');
title('Fault-tolerant mode activation (Scenario 3)');

%% ================== Report metrics (steady window) ==================
win = time > (Tend-0.5);
report_metrics('Healthy', res1, win);
report_metrics('Fault w/o FTC', res2, win);
report_metrics('Fault WITH FTC', res3, win);

% speed estimation quality
fprintf('\n=== Speed estimation error (last 0.5 s) ===\n');
print_speed_err('Healthy', res1, win);
print_speed_err('Fault w/o FTC', res2, win);
print_speed_err('Fault WITH FTC', res3, win);

%% ================== Helpers ==================
function opt = make_opts(name, enable_fault, enable_ftc)
    opt.name = name;
    opt.enable_fault = enable_fault;
    opt.enable_ftc = enable_ftc;
end

function report_metrics(tag, res, win)
    Te_mean = mean(res.Te(win));
    Te_pp   = max(res.Te(win)) - min(res.Te(win));
    w_mean  = mean(res.wm_rpm(win));
    w_pp    = max(res.wm_rpm(win)) - min(res.wm_rpm(win));
    ia_r = rms(res.ia(win));
    ib_r = rms(res.ib(win));
    ic_r = rms(res.ic(win));
    ratio_ia = ia_r / max(ib_r,1e-12);
    ratio_sum = rms(res.ib(win)+res.ic(win)) / max(ib_r,1e-12);

    fprintf('\n%s:\n', tag);
    fprintf('Te_mean=%.2f N.m, Te_pp=%.2f N.m\n', Te_mean, Te_pp);
    fprintf('w_mean=%.2f rpm, w_ripple_pp=%.2f rpm\n', w_mean, w_pp);
    fprintf('Ia_rms/Ib_rms=%.4f\n', ratio_ia);
    fprintf('RMS(ib+ic)/Ib_rms=%.4f\n', ratio_sum);
end

function print_speed_err(tag, res, win)
    e = res.wm_hat_rpm(win) - res.wm_rpm(win);
    fprintf('%s: mean|err|=%.2f rpm, max|err|=%.2f rpm\n', tag, mean(abs(e)), max(abs(e)));
end

%% ================== Main simulation ==================
function out = sim_case(opt)
    % pull from base
    Rs = evalin('base','Rs');
    Rr = evalin('base','Rr');
    Lm = evalin('base','Lm');
    Ls = evalin('base','Ls');
    Lr = evalin('base','Lr');
    Delta = evalin('base','Delta');
    Tr = evalin('base','Tr');
    P = evalin('base','P');
    J = evalin('base','J');
    B = evalin('base','B');

    Ts = evalin('base','Ts');
    time = evalin('base','time');
    N = evalin('base','N');
    tfault = evalin('base','tfault');

    TL = evalin('base','TL');
    w_ref = evalin('base','w_ref');
    ws_ref = evalin('base','ws_ref');

    Vpk_max = evalin('base','Vpk_max');
    Kvf = evalin('base','Kvf');
    Kvf_fault = evalin('base','Kvf_fault');

    Kpw = evalin('base','Kpw');
    Kiw = evalin('base','Kiw');
    Kpw_fault = evalin('base','Kpw_fault');
    Kiw_fault = evalin('base','Kiw_fault');

    ws_min = evalin('base','ws_min');
    ws_max = evalin('base','ws_max');

    % MRAS
    kf_flux = evalin('base','kf_flux');
    Kp_mras = evalin('base','Kp_mras');
    Ki_mras = evalin('base','Ki_mras');
    wr_hat_min = evalin('base','wr_hat_min');
    wr_hat_max = evalin('base','wr_hat_max');

    % FDI
    Nw = evalin('base','Nw');
    thr_ia  = evalin('base','thr_ia');
    thr_sum = evalin('base','thr_sum');
    hold_windows = evalin('base','hold_windows');

    % plant states
    psi_s = 0 + 1j*0;
    psi_r = 0 + 1j*0;
    wm = 0;

    % controller states
    theta = 0;
    int_w = 0;

    % sensorless MRAS states
    psi_s_v = 0 + 1j*0;   % voltage-model stator flux
    psi_r_i = 0 + 1j*0;   % current-model rotor flux
    int_mras = 0;
    wr_hat = 0;           % electrical rotor speed estimate

    % FDI buffers (running RMS)
    ia_buf = zeros(1,Nw);
    ib_buf = zeros(1,Nw);
    ic_buf = zeros(1,Nw);
    sum_buf = zeros(1,Nw);
    sum_ia2 = 0; sum_ib2 = 0; sum_sum2 = 0;
    idx = 1;
    fcnt = 0;
    ftc_active = false;

    % logs
    ia = zeros(1,N); ib = zeros(1,N); ic = zeros(1,N);
    wm_rpm = zeros(1,N);
    wm_hat_rpm = zeros(1,N);
    Te_log = zeros(1,N);
    ftc_log = zeros(1,N);

    for k = 1:N
        t = time(k);
        omega_r = (P/2)*wm;

        fault = opt.enable_fault && (t >= tfault);

        % ===== Enforce open-phase constraint at start of step =====
        if fault
            psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
        end

        % ===== Currents from plant state =====
        i_s = (Lr*psi_s - Lm*psi_r)/Delta;
        i_alpha = real(i_s);
        i_beta  = imag(i_s);

        % Phase currents for measurement/plots
        if fault
            ia_k = 0;
            ib_k =  (sqrt(3)/2)*i_beta;
            ic_k = -(sqrt(3)/2)*i_beta;
        else
            ia_k = i_alpha;
            ib_k = -0.5*i_alpha + (sqrt(3)/2)*i_beta;
            ic_k = -0.5*i_alpha - (sqrt(3)/2)*i_beta;
        end
        ia(k)=ia_k; ib(k)=ib_k; ic(k)=ic_k;

        % ===== Fault Detection (FDI) =====
        % Update running sums (ring buffer)
        ia_old = ia_buf(idx); ib_old = ib_buf(idx); ic_old = ic_buf(idx); sum_old = sum_buf(idx);
        sum_ia2 = sum_ia2 - ia_old^2;  sum_ib2 = sum_ib2 - ib_old^2;  sum_sum2 = sum_sum2 - sum_old^2;

        ia_buf(idx) = ia_k; ib_buf(idx) = ib_k; ic_buf(idx) = ic_k;
        s_k = ib_k + ic_k;
        sum_buf(idx) = s_k;

        sum_ia2 = sum_ia2 + ia_k^2;
        sum_ib2 = sum_ib2 + ib_k^2;
        sum_sum2 = sum_sum2 + s_k^2;

        idx = idx + 1; if idx > Nw, idx = 1; end

        if k >= Nw
            ia_rms = sqrt(sum_ia2/Nw);
            ib_rms = sqrt(sum_ib2/Nw);
            sum_rms = sqrt(sum_sum2/Nw);

            r_ia  = ia_rms / max(ib_rms,1e-12);
            r_sum = sum_rms / max(ib_rms,1e-12);

            cond = (r_ia < thr_ia) && (r_sum < thr_sum) && (t > 0.5); % avoid early startup
            if cond
                fcnt = fcnt + 1;
            else
                fcnt = 0;
            end

            if opt.enable_ftc && ~ftc_active && (fcnt >= hold_windows)
                ftc_active = true;
                % Reconfiguration actions (typical supervisory FTC)
                int_w = 0;           % reset speed PI integrator
                int_mras = 0;        % reset MRAS integrator
            end
        end

        ftc_log(k) = ftc_active;

        % ===== Sensorless MRAS estimator =====
        % We use applied v_s (computed below), so compute control first with last wr_hat.

        % Estimated mechanical speed for control (use wr_hat)
        wm_hat = (2/P)*wr_hat;  % mechanical rad/s

        % ===== Speed PI (uses estimated speed) =====
        if ftc_active
            Kpw_use = Kpw_fault; Kiw_use = Kiw_fault; Kvf_use = Kvf_fault;
        else
            Kpw_use = Kpw; Kiw_use = Kiw; Kvf_use = Kvf;
        end

        e_w = w_ref - wm_hat;
        int_w = int_w + e_w*Ts;
        wslip = Kpw_use*e_w + Kiw_use*int_w;
        ws_unsat = ws_ref + wslip;
        ws_cmd = min(max(ws_unsat, ws_min), ws_max);
        if ws_cmd ~= ws_unsat
            int_w = int_w + (ws_cmd - ws_unsat)/max(Kiw_use,1e-9);
        end

        theta = theta + ws_cmd*Ts;

        % ===== Voltage synthesis =====
        Vpk = Kvf_use * ws_cmd;
        Vpk = min(max(Vpk, 0), Vpk_max);

        if ~fault
            % Healthy: balanced phase voltages
            va_cmd = Vpk*sin(theta);
            vb_cmd = Vpk*sin(theta - 2*pi/3);
            vc_cmd = Vpk*sin(theta + 2*pi/3);

            vn = (va_cmd + vb_cmd + vc_cmd)/3;
            va = va_cmd - vn;
            vb = vb_cmd - vn;
            vc = vc_cmd - vn;

            v_alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
            v_beta  = (2/3)*((sqrt(3)/2)*(vb - vc));
        else
            % Whole-leg open (phase A leg). Plant is driven by line voltage v_bc.
            if ftc_active
                % FTC mode: direct two-leg modulation (use a clean sinus for v_bc)
                % Choose v_bc so that resulting v_beta amplitude equals Vpk
                v_bc = sqrt(3)*Vpk*sin(theta);
            else
                % No FTC: keep using pre-fault balanced vb,vc commands (va leg is unavailable)
                vb_cmd = Vpk*sin(theta - 2*pi/3);
                vc_cmd = Vpk*sin(theta + 2*pi/3);
                v_bc = vb_cmd - vc_cmd;  % equals -sqrt(3)*Vpk*cos(theta)
            end

            v_alpha = 0;
            v_beta  = (sqrt(3)/3)*v_bc;
        end

        v_s = v_alpha + 1j*v_beta;

        % ===== MRAS update (needs v_s and i_s) =====
        % voltage model stator flux with leakage correction
        psi_s_v = psi_s_v + Ts*(v_s - Rs*i_s - kf_flux*psi_s_v);

        % current model rotor flux
        dpsi_r_i = (Lm/Tr)*i_s - (1/Tr)*psi_r_i + 1j*wr_hat*psi_r_i;
        psi_r_i = psi_r_i + Ts*dpsi_r_i;

        % current model stator flux
        psi_s_i = (Delta/Lr)*i_s + (Lm/Lr)*psi_r_i;

        % normalized MRAS error (flux angle mismatch)
        denom = max(abs(psi_s_v)^2, 1e-6);
        e_mras = imag(conj(psi_s_i)*psi_s_v) / denom;

        int_mras = int_mras + e_mras*Ts;
        wr_hat = wr_hat + Ts*(Kp_mras*e_mras + Ki_mras*int_mras);
        wr_hat = min(max(wr_hat, wr_hat_min), wr_hat_max);

        wm_hat_rpm(k) = (2/P)*wr_hat * 60/(2*pi);

        % ===== Plant flux dynamics =====
        dpsi_s = v_s - Rs*i_s;
        i_r = (-Lm*psi_s + Ls*psi_r)/Delta;
        dpsi_r = -Rr*i_r + 1j*omega_r*psi_r;

        psi_s = psi_s + dpsi_s*Ts;
        psi_r = psi_r + dpsi_r*Ts;

        % re-apply open-phase constraint
        if fault
            psi_s = (Lm/Lr)*real(psi_r) + 1j*imag(psi_s);
            i_s = (Lr*psi_s - Lm*psi_r)/Delta;
        end

        % torque + mech
        Te = 1.5*(P/2) * imag(conj(psi_s) * i_s);
        wm = wm + (Te - TL - B*wm)/J * Ts;

        wm_rpm(k) = wm*60/(2*pi);
        Te_log(k) = Te;
    end

    out.name = opt.name;
    out.wm_rpm = wm_rpm;
    out.wm_hat_rpm = wm_hat_rpm;
    out.Te = Te_log;
    out.ia = ia; out.ib = ib; out.ic = ic;
    out.ftc_active = ftc_log;
end
