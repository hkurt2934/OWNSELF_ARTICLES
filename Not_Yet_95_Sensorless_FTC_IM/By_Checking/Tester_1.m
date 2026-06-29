%% =========================================================
% DIRECT FOC + MRAS SPEED ESTIMATOR
% Induction Motor - SVPWM Average Inverter Model
%% =========================================================

clear; clc;

%% ===============================
% 1. Motor Parameters
%% ===============================
Rs = 1.405;
Rr = 1.395;
Ls = 0.0058;
Lr = 0.0058;
Lm = 0.0055;

P  = 2;        % pole pairs
J  = 0.02;
B  = 0.001;

Tr = Lr/Rr;

%% ===============================
% 2. Simulation Parameters
%% ===============================
Ts   = 1e-4;
Tsim = 3;
N    = round(Tsim/Ts);
time = (0:N-1)*Ts;

Vdc  = 560;
Vmax = Vdc/sqrt(3);    % SVPWM limit
TL   = 10;

%% ===============================
% 3. Control Parameters
%% ===============================
% Speed controller
Kpw = 1.5;
Kiw = 50;

% Current controllers
Kpd = 40;   Kid = 400;
Kpq = 40;   Kiq = 400;

% Anti-windup
Kaw_d = 1/Kpd;
Kaw_q = 1/Kpq;

% MRAS
Kp_mras = 200;
Ki_mras = 4000;

%% ===============================
% 4. Initialization
%% ===============================
ids = 0; iqs = 0;
wr  = 0;

theta_e = 0;

% Integrators
int_w  = 0;
int_id = 0;
int_iq = 0;
int_mras = 0;

% Stator flux observer (Direct FOC)
psi_s_alpha = 0;
psi_s_beta  = 0;

% MRAS
psi_ref = [0;0];
psi_ad  = [0;0];
w_e = 0;

%% ===============================
% References
%% ===============================
w_ref = 100;           % rad/s
psi_s_ref = 0.8;
ids_ref = psi_s_ref/Lm;

%% ===============================
% Logging
%% ===============================
wr_log = zeros(1,N);
wr_hat_log = zeros(1,N);
Te_log = zeros(1,N);

ids_log = zeros(1,N);
iqs_log = zeros(1,N);
psi_s_mag_log = zeros(1,N);
speed_err_log = zeros(1,N);

%% ===============================
% Main Simulation Loop
%% ===============================
for k = 1:N

    %% -------- Speed Controller --------
    ew = w_ref - w_e;
    int_w = int_w + ew*Ts;
    iqs_ref = Kpw*ew + Kiw*int_w;

    %% -------- Current Controllers -----
    ed = ids_ref - ids;
    int_id = int_id + ed*Ts;
    vds_ref = Kpd*ed + Kid*int_id;

    eq = iqs_ref - iqs;
    int_iq = int_iq + eq*Ts;
    vqs_ref = Kpq*eq + Kiq*int_iq;

    %% -------- SVPWM DC-Link Constraint --------
    Vref = sqrt(vds_ref^2 + vqs_ref^2);

    if Vref > Vmax
        vds = vds_ref * Vmax/Vref;
        vqs = vqs_ref * Vmax/Vref;
    else
        vds = vds_ref;
        vqs = vqs_ref;
    end

    % Anti-windup
    int_id = int_id + Kaw_d*(vds - vds_ref)*Ts;
    int_iq = int_iq + Kaw_q*(vqs - vqs_ref)*Ts;

    %% -------- Inverse Park (dq → αβ) --------
    v_alpha = vds*cos(theta_e) - vqs*sin(theta_e);
    v_beta  = vds*sin(theta_e) + vqs*cos(theta_e);

    %% -------- Electromagnetic Torque --------
    Te = (3/2)*P*(Lm/Lr)*psi_s_ref*iqs;

    %% -------- Mechanical Dynamics --------
    dwr = (Te - TL - B*wr)/J;
    wr = wr + dwr*Ts;

    %% -------- Stator Current Dynamics --------
    dids = (vds - Rs*ids + w_e*Ls*iqs)/Ls;
    diqs = (vqs - Rs*iqs - w_e*Ls*ids)/Ls;

    ids = ids + Ts*dids;
    iqs = iqs + Ts*diqs;

    %% -------- Stator Flux Observer (Direct FOC) --------
    i_alpha = ids*cos(theta_e) - iqs*sin(theta_e);
    i_beta  = ids*sin(theta_e) + iqs*cos(theta_e);

    psi_s_alpha = psi_s_alpha + Ts*(v_alpha - Rs*i_alpha);
    psi_s_beta  = psi_s_beta  + Ts*(v_beta  - Rs*i_beta);

    theta_e = atan2(psi_s_beta, psi_s_alpha);

    %% -------- MRAS Speed Estimator --------
    % Voltage model
    psi_ref(1) = psi_ref(1) + Ts*(v_alpha - Rs*i_alpha);
    psi_ref(2) = psi_ref(2) + Ts*(v_beta  - Rs*i_beta);

    % Current model
    psi_ad(1) = psi_ad(1) + Ts*(-psi_ad(1)/Tr + w_e*psi_ad(2) + (Lm/Tr)*i_alpha);
    psi_ad(2) = psi_ad(2) + Ts*(-psi_ad(2)/Tr - w_e*psi_ad(1) + (Lm/Tr)*i_beta);

    e_mras = psi_ref(1)*psi_ad(2) - psi_ref(2)*psi_ad(1);
    int_mras = int_mras + e_mras*Ts;
    w_e = Kp_mras*e_mras + Ki_mras*int_mras;

    %% -------- Logging --------
    wr_log(k) = wr;
    wr_hat_log(k) = w_e;
    Te_log(k) = Te;

    ids_log(k) = ids;
    iqs_log(k) = iqs;
    psi_s_mag_log(k) = sqrt(psi_s_alpha^2 + psi_s_beta^2);
    speed_err_log(k) = w_ref - w_e;

end

%% ===============================
% PLOTTING RESULTS
%% ===============================

% Speed
subplot(3,2,1);
plot(time, wr_log, 'LineWidth',1.5); hold on;
plot(time, wr_hat_log,'--','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Speed (rad/s)');
legend('Actual Speed','Estimated Speed (MRAS)');
grid on;
title('Speed Response');

% Speed error
subplot(3,2,2);
plot(time, speed_err_log,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Speed Error (rad/s)');
grid on;
title('Speed Estimation Error');

% Torque
subplot(3,2,3);
plot(time, Te_log,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Torque (Nm)');
grid on;
title('Electromagnetic Torque');

% dq currents
subplot(3,2,4);
plot(time, ids_log,'LineWidth',1.5); hold on;
plot(time, iqs_log,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Current (A)');
legend('i_d','i_q');
grid on;
title('dq-Axis Currents');

% Stator flux magnitude
subplot(3,2,5);
plot(time, psi_s_mag_log,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('|\psi_s| (Wb)');
grid on;
title('Stator Flux Magnitude');
