%% ============================================================
%  SENSORLESS INDUCTION MOTOR DRIVE WITH FOC AND FTC
% ============================================================
clear; clc; close all;

%% ---------------- SIMULATION PARAMETERS ----------------
Ts   = 1e-4;          % Sampling time [s]
Tsim = 3.0;           % Simulation time [s]
N    = round(Tsim/Ts);
t    = (0:N-1)*Ts;

%% ---------------- MOTOR PARAMETERS ----------------
Rs = 1.405;           % Stator resistance [Ohm]
Rr = 1.395;           % Rotor resistance [Ohm]
Ls = 0.0058;          % Stator inductance [H]
Lr = 0.0058;          % Rotor inductance [H]
Lm = 0.0055;          % Magnetizing inductance [H]
P  = 2;               % Pole pairs
J  = 0.02;            % Inertia [kg.m^2]
B  = 0.001;           % Friction coefficient

sigma = 1 - (Lm^2)/(Ls*Lr);

%% ---------------- DC LINK & LIMITS ----------------
Vdc  = 300;
Vmax = Vdc/sqrt(3);
Imax = 10;

%% ---------------- REFERENCES ----------------
w_ref = 150;           % Speed reference [rad/s]
id_ref = 2;            % Flux current reference [A]

%% ---------------- PI CONTROLLER GAINS ----------------
% Speed loop (slow)
Kp_w = 1.0;
Ki_w = 50;

% Current loop (fast)
Kp_i = 20;
Ki_i = 2000;

%% ---------------- INITIAL CONDITIONS ----------------
id = 0; iq = 0;
w_m = 0;
theta_e = 0;

psi_r = 0.8;           % Rotor flux initial value

int_w  = 0;
int_id = 0;
int_iq = 0;

%% ---------------- FAULT SETUP ----------------
t_fault = 1.0;         % Fault time
fault_gain = 1;

%% ---------------- DATA STORAGE ----------------
w_log    = zeros(1,N);
we_log   = zeros(1,N);
id_log   = zeros(1,N);
iq_log   = zeros(1,N);
Te_log   = zeros(1,N);
psi_log  = zeros(1,N);
err_log  = zeros(1,N);

%% ====================================================
%                    MAIN LOOP
% ====================================================
for k = 1:N

    %% -------- FAULT INJECTION (Rs INCREASE) --------
    if t(k) >= t_fault
        fault_gain = 1.5;     % 50% stator resistance fault
    end
    Rs_eff = Rs * fault_gain;

    %% -------- SPEED CONTROL LOOP --------
    w_err = w_ref - w_m;
    int_w = int_w + Ki_w * Ts * w_err;
    int_w = max(min(int_w, Imax), -Imax);

    iq_ref = Kp_w * w_err + int_w;
    iq_ref = max(min(iq_ref, Imax), -Imax);

    %% -------- CURRENT CONTROL LOOPS --------
    err_id = id_ref - id;
    err_iq = iq_ref - iq;

    int_id = int_id + Ki_i * Ts * err_id;
    int_iq = int_iq + Ki_i * Ts * err_iq;

    vd = Kp_i * err_id + int_id;
    vq = Kp_i * err_iq + int_iq;

    %% -------- VOLTAGE LIMITATION --------
    Vmag = sqrt(vd^2 + vq^2);
    if Vmag > Vmax
        vd = vd * Vmax / Vmag;
        vq = vq * Vmax / Vmag;
    end

    %% -------- ELECTRICAL DYNAMICS --------
    did = (vd - Rs_eff*id + w_m*Ls*iq) / Ls;
    diq = (vq - Rs_eff*iq - w_m*Ls*id) / Ls;

    id = id + Ts * did;
    iq = iq + Ts * diq;

    %% -------- ROTOR FLUX OBSERVER --------
    dpsi = (Lm/Lr)*id - (Rr/Lr)*psi_r;
    psi_r = psi_r + Ts * dpsi;
    psi_r = max(psi_r, 0.2);

    %% -------- ELECTROMAGNETIC TORQUE --------
    Te = (3/2) * P * (Lm/Lr) * psi_r * iq;

    %% -------- MECHANICAL DYNAMICS --------
    
    dw = (Te - B*w_m) / J;
    w_m = w_m + Ts * dw;

    %% -------- ELECTRICAL ANGLE --------
    theta_e = theta_e + Ts * P * w_m;

    %% -------- LOGGING --------
    w_log(k)   = w_m;
    id_log(k)  = id;
    iq_log(k)  = iq;
    Te_log(k)  = Te;
    psi_log(k) = psi_r;
    err_log(k) = w_err;

end

%% ====================================================
%                       PLOTS
% ====================================================

% Rotor Speed
figure;
plot(t, w_log, 'b', 'LineWidth',1.8); hold on;
yline(w_ref,'k--','LineWidth',1.5);
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('Speed (rad/s)');
title('Rotor Speed Response');
grid on;

% Speed Error
figure;
plot(t, err_log, 'r', 'LineWidth',1.8);
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('Speed Error (rad/s)');
title('Speed Tracking Error');
grid on;

% iq Current
figure;
plot(t, iq_log, 'b', 'LineWidth',1.8); hold on;
yline(mean(iq_log(1:200)),'k--');
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('i_q (A)');
title('q-axis Current');
grid on;

% id Current
figure;
plot(t, id_log, 'b', 'LineWidth',1.8); hold on;
yline(id_ref,'k--','i_d^*');
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('i_d (A)');
title('d-axis Current');
grid on;

% Electromagnetic Torque
figure;
plot(t, Te_log, 'b', 'LineWidth',1.8);
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('Torque (Nm)');
title('Electromagnetic Torque');
grid on;

% Rotor Flux
figure;
plot(t, psi_log, 'b', 'LineWidth',1.8);
xline(t_fault,'--','Fault');
xlabel('Time (s)');
ylabel('\psi_r (Wb)');
title('Rotor Flux Magnitude');
grid on;
