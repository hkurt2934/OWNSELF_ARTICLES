%% Step 2 - Healthy IM + FOC with SPWM (average VSI)
clear; clc; close all;

%% Timing
Ts = 100e-6; Tsim = 3.0; N = round(Tsim/Ts)+1; t = (0:N-1)'*Ts;
t_ramp = 0.3;                           % [s]

%% Motor parameters (typical 1.5 kW, 4-pole, 400 V class)
P = 4;  p = P/2;
Rs = 2.2;  Rr = 1.6;
Ls = 0.21; Lr = 0.21; Lm = 0.20;
J  = 0.02; B  = 0.0;

% Derived
sigma = 1 - (Lm^2)/(Ls*Lr);
Lsig  = sigma*Ls;
Tr = Lr/Rr;
Kt = 1.5*p*(Lm/Lr); % Te = Kt * psi_r * iqs (psi_r aligned)

%% DC link + limits
Vdc = 560;               % you used this earlier (reasonable)
Vph_max = 0.48*Vdc;       % SPWM linear-ish phase peak limit (~Vdc/2)
Imax = 20;                % current limit

%% References
w_ref_rpm = 1500; 
w_ref = w_ref_rpm*2*pi/60;
Tload = 9.55;             % if you want rated: ~9.55 N.m for 1.5kW@1500rpm
psi_r_ref = 0.85;         % Wb (flux reference)
ids_ref = psi_r_ref/Lm;   % A (approx in aligned frame)
ids_ref = min(max(ids_ref, 0), Imax);

%% Controllers
% Current loop bandwidth
wc_i = 2*pi*300;          % rad/s
Kp_id = Lsig*wc_i; Ki_id = Rs*wc_i;
Kp_iq = Lsig*wc_i; Ki_iq = Rs*wc_i;

% Speed loop (slow)
wc_w = 2*pi*15;
Kp_w = 1;               % tuned conservative
Ki_w = 12;

% Anti-windup gains (back-calculation)
Kaw_i = 200;  Kaw_w = 20;

%% State init (plant currents in dq synchronous frame)
ids=0; iqs=0; idr=0; iqr=0;
wm=0;                     % mech rad/s

% Rotor flux estimate for control (scalar along d-axis)
psi_r_hat = 0.01;

% Angles
theta_e = 0;

% Integrators
xi_w=0; xi_id=0; xi_iq=0;

%% Logs
wm_log=zeros(N,1); Te_log=zeros(N,1);
ids_log=zeros(N,1); iqs_log=zeros(N,1);
idsr_log=zeros(N,1); iqsr_log=zeros(N,1);
we_log=zeros(N,1); wsl_log=zeros(N,1);
da_log=zeros(N,1); db_log=zeros(N,1); dc_log=zeros(N,1);

%% Helper functions
wrap = @(x) mod(x+pi,2*pi)-pi;

for k=1:N
    % Electrical rotor speed
    wr_e = p*wm;

    % --- Rotor flux estimator (current model, aligned) ---
    % psi_r_dot = (Lm/Tr)*ids - (1/Tr)*psi_r
    psi_r_hat = psi_r_hat + Ts*((Lm/Tr)*ids - (1/Tr)*psi_r_hat);
    if abs(psi_r_hat) < 1e-3, psi_r_hat = sign(psi_r_hat+1e-6)*1e-3; end

    % --- Speed PI -> torque current reference ---
    ew = w_ref - wm;
    Te_cmd_unsat = Kp_w*ew + xi_w;
    % Convert torque to iqs* using Te = Kt*psi_r*i_qs
    iqs_ref_unsat = Te_cmd_unsat/(Kt*psi_r_hat);
    iqs_ref = min(max(iqs_ref_unsat, -Imax), Imax);

    % Speed PI anti-windup
    xi_w = xi_w + Ts*(Ki_w*ew + Kaw_w*(iqs_ref - iqs_ref_unsat));

    % --- Slip + synchronous speed ---
    w_sl = (Lm/Tr)*(iqs_ref/psi_r_hat);
    we = wr_e + w_sl;
    theta_e = wrap(theta_e + Ts*we);

    % --- Current control ---
    eid = ids_ref - ids;
    eiq = iqs_ref - iqs;

    vds_u = Kp_id*eid + xi_id;
    vqs_u = Kp_iq*eiq + xi_iq;

    % Feedforward / decoupling (simple, helps stability)
    vds_ff = Rs*ids - we*Lsig*iqs;
    vqs_ff = Rs*iqs + we*(Lsig*ids + (Lm/Lr)*psi_r_hat);

    vds_star_unsat = vds_u + vds_ff;
    vqs_star_unsat = vqs_u + vqs_ff;

    % Voltage vector saturation in alpha-beta (equivalent)
    % Convert to alpha-beta using inverse Park with theta_e
    v_alpha =  cos(theta_e)*vds_star_unsat - sin(theta_e)*vqs_star_unsat;
    v_beta  =  sin(theta_e)*vds_star_unsat + cos(theta_e)*vqs_star_unsat;
    Vmag = hypot(v_alpha, v_beta);
    if Vmag > Vph_max
        scale = Vph_max/Vmag;
        v_alpha = v_alpha*scale;
        v_beta  = v_beta*scale;
    end
    % Back to dq after saturation
    vds_star =  cos(theta_e)*v_alpha + sin(theta_e)*v_beta;
    vqs_star = -sin(theta_e)*v_alpha + cos(theta_e)*v_beta;

    % Current PI anti-windup (back-calc with actual sat effect)
    xi_id = xi_id + Ts*(Ki_id*eid + Kaw_i*(vds_star - vds_star_unsat));
    xi_iq = xi_iq + Ts*(Ki_iq*eiq + Kaw_i*(vqs_star - vqs_star_unsat));

    % --- SPWM average inverter ---
    % alpha-beta -> abc phase commands
    va = v_alpha;
    vb = -0.5*v_alpha + (sqrt(3)/2)*v_beta;
    vc = -0.5*v_alpha - (sqrt(3)/2)*v_beta;

    % duties: v_phase = (d-0.5)*Vdc
    da = 0.5 + va/Vdc;
    db = 0.5 + vb/Vdc;
    dc = 0.5 + vc/Vdc;
    da = min(max(da,0),1);
    db = min(max(db,0),1);
    dc = min(max(dc,0),1);

    % Applied phase voltages (average)
    va_ap = (da-0.5)*Vdc;
    vb_ap = (db-0.5)*Vdc;
    vc_ap = (dc-0.5)*Vdc;

    % Applied alpha-beta
    v_alpha_ap = (2/3)*(va_ap - 0.5*vb_ap - 0.5*vc_ap);
    v_beta_ap  = (2/3)*((sqrt(3)/2)*(vb_ap - vc_ap));

    % Applied dq
    vds =  cos(theta_e)*v_alpha_ap + sin(theta_e)*v_beta_ap;
    vqs = -sin(theta_e)*v_alpha_ap + cos(theta_e)*v_beta_ap;

    % --- Plant model (dq currents + rotor currents) ---
    % Using current state model (same structure as many dq IM models)
    % x = [ids iqs idr iqr]
    % L matrix
    % [vds]   [Rs 0  0  0][ids] + d/dt([Ls 0 Lm 0][ids]) + cross terms
    % For numerical robustness: solve flux-based derivatives

    % Fluxes
    psi_ds = Ls*ids + Lm*idr;
    psi_qs = Ls*iqs + Lm*iqr;
    psi_dr = Lm*ids + Lr*idr;
    psi_qr = Lm*iqs + Lr*iqr;

    % Flux derivatives (Krause in synchronous frame)
    dpsi_ds = vds - Rs*ids + we*psi_qs;
    dpsi_qs = vqs - Rs*iqs - we*psi_ds;

    dpsi_dr = -Rr*idr + (we - wr_e)*psi_qr;
    dpsi_qr = -Rr*iqr - (we - wr_e)*psi_dr;

    % Convert dpsi -> di using inverse inductance matrix
    % [psi_s] = [Ls Lm; Lm Lr] [i_s; i_r]
    % Invert 2x2 for each axis (d and q)
    detL = Ls*Lr - Lm*Lm;
    ids_dot = ( Lr*dpsi_ds - Lm*dpsi_dr)/detL;
    idr_dot = (-Lm*dpsi_ds + Ls*dpsi_dr)/detL;
    iqs_dot = ( Lr*dpsi_qs - Lm*dpsi_qr)/detL;
    iqr_dot = (-Lm*dpsi_qs + Ls*dpsi_qr)/detL;

    ids = ids + Ts*ids_dot;
    iqs = iqs + Ts*iqs_dot;
    idr = idr + Ts*idr_dot;
    iqr = iqr + Ts*iqr_dot;

    % Electromagnetic torque
    Te = 1.5*p*(psi_ds*iqs - psi_qs*ids);

    % Mechanical
    wm = wm + Ts*((Te - Tload - B*wm)/J);

    % Log
    wm_log(k)=wm;
    Te_log(k)=Te;
    ids_log(k)=ids; iqs_log(k)=iqs;
    idsr_log(k)=ids_ref; iqsr_log(k)=iqs_ref;
    we_log(k)=we; wsl_log(k)=w_sl;
    da_log(k)=da; db_log(k)=db; dc_log(k)=dc;
end

%% Plots
wm_rpm = wm_log*60/(2*pi);

figure('Name','Speed');
plot(t, w_ref_rpm*ones(size(t)),'LineWidth',1.2); hold on;
plot(t, wm_rpm,'LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('Speed [rpm]');
legend('Reference','Actual'); title('FOC Speed Control (Healthy, SPWM avg)');

figure('Name','dq currents');
plot(t, ids_log,'LineWidth',1.1); hold on;
plot(t, iqs_log,'LineWidth',1.1);
plot(t, idsr_log,'--','LineWidth',1.1);
plot(t, iqsr_log,'--','LineWidth',1.1);
grid on; xlabel('Time [s]'); ylabel('Current [A]');
legend('i_{ds}','i_{qs}','i_{ds}^*','i_{qs}^*'); title('Stator dq currents');

figure('Name','Torque');
plot(t, Te_log,'LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('T_e [N*m]'); title('Electromagnetic torque');

figure('Name','SPWM duty (first 40 ms)');
idx = t<=0.04;
plot(t(idx), da_log(idx),'LineWidth',1.1); hold on;
plot(t(idx), db_log(idx),'LineWidth',1.1);
plot(t(idx), dc_log(idx),'LineWidth',1.1);
grid on; xlabel('Time [s]'); ylabel('Duty'); legend('d_a','d_b','d_c');
