clc;
clear;
close all;

%% =========================
% Simulation Parameters
% =========================
Ts   = 1e-4;        % Sampling time [s]
Tsim = 3.0;         % Simulation time [s]
t = 0:Ts:Tsim;

%% =========================
% Induction Motor Parameters
% =========================
Rs = 1.405;         % Stator resistance [ohm]
Rr = 1.395;         % Rotor resistance [ohm]
Ls = 0.0058;        % Stator inductance [H]
Lr = 0.0058;        % Rotor inductance [H]
Lm = 0.0055;        % Mutual inductance [H]
P  = 2;             % Pole pairs

J = 0.01;           % Inertia [kg.m^2]
B = 0.001;          % Friction coefficient

%% =========================
% Inverter / DC Link
% =========================
Vdc = 560;          % DC link voltage [V]

%% =========================
% Load Torque
% =========================
TL = 10;            % Constant load torque [N.m]

%% =========================
% Derived Parameters
% =========================
sigma = 1 - (Lm^2)/(Ls*Lr);
Tr = Lr / Rr;

%% =========================
% Control Gains
% =========================
Kp_id = 20;     Ki_id = 2000;
Kp_iq = 20;     Ki_iq = 2000;
Kp_w  = 2.5;    Ki_w  = 100;

%% =========================
% References
% =========================
w_ref = 1200 * 2*pi/60;     % Speed reference [rad/s]
psi_r_ref = 0.8;            % Rotor flux reference [Wb]

%% =========================
% Initial Conditions
% =========================
id = 0;  iq = 0;
psi_r = 0.1;

omega_m = 0;
theta_e = 0;

int_id = 0;
int_iq = 0;
int_w  = 0;

%% =========================
% Data Logging
% =========================
id_log = zeros(size(t));
iq_log = zeros(size(t));
w_log  = zeros(size(t));
Te_log = zeros(size(t));

%% =========================
% Main Simulation Loop
% =========================
for k = 1:length(t)

    %% Speed Controller
    
    w_err = w_ref - omega_m;
    int_w = int_w + w_err * Ts;
    iq_ref = Kp_w * w_err + Ki_w * int_w;

    %% Flux Controller (d-axis current reference)
    id_ref = psi_r_ref / Lm;

    %% d-axis Current Controller
    id_err = id_ref - id;
    int_id = int_id + id_err * Ts;
    vd = Kp_id * id_err + Ki_id * int_id;

    %% q-axis Current Controller
    iq_err = iq_ref - iq;
    int_iq = int_iq + iq_err * Ts;
    vq = Kp_iq * iq_err + Ki_iq * int_iq;

    %% Voltage Limitation (SPWM linear region)
    Vmax = Vdc / sqrt(3);
    Vmag = sqrt(vd^2 + vq^2);
    if Vmag > Vmax
        vd = vd * Vmax / Vmag;
        vq = vq * Vmax / Vmag;
    end

    %% Rotor Flux Dynamics
    dpsi_r = (Lm * id - psi_r) / Tr;
    psi_r = psi_r + dpsi_r * Ts;

    %% Slip & Electrical Speed
    omega_sl = (Lm * Rr / Lr) * (iq / psi_r);
    omega_e  = P * omega_m + omega_sl;

    %% Induction Motor dq Model
    did = (vd - Rs*id + omega_e*sigma*Ls*iq) / (sigma*Ls);
    diq = (vq - Rs*iq - omega_e*sigma*Ls*id ...
          - omega_e*(Lm/Lr)*psi_r) / (sigma*Ls);

    id = id + did * Ts;
    iq = iq + diq * Ts;

    %% Electromagnetic Torque
    Te = (3/2) * P * (Lm/Lr) * psi_r * iq;

    %% Mechanical Dynamics (with load torque)
    domega = (Te - TL - B*omega_m) / J;
    omega_m = omega_m + domega * Ts;

    %% Electrical Angle
    theta_e = theta_e + omega_e * Ts;
    theta_e = mod(theta_e, 2*pi);

    %% Logging
    id_log(k) = id;
    iq_log(k) = iq;
    w_log(k)  = omega_m;
    Te_log(k) = Te;

end

%% =========================
% Plot Results
% =========================
figure;
plot(t, w_log*60/(2*pi), 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Speed [rpm]');
title('Induction Motor Speed Response (FOC, TL = 10 N.m)');

figure;
plot(t, id_log, 'LineWidth', 1.6); hold on;
plot(t, iq_log, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Current [A]');
legend('i_d','i_q');
title('dq-axis Currents');

figure;
plot(t, Te_log, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Torque [N.m]');
title('Electromagnetic Torque');
