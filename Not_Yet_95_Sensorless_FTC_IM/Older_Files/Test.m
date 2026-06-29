%% SENSORLESS INDUCTION MOTOR DRIVE WITH FTC
clear; clc; close all;

%% ---------------- SIMULATION PARAMETERS ----------------
Ts   = 1e-4;              % Sampling time [s]
Tsim = 3.0;               % Total simulation time [s]
N    = round(Tsim/Ts);
t    = (0:N-1)*Ts;

%% ---------------- MOTOR PARAMETERS ----------------
Rs = 1.405;
Rr = 1.395;
Ls = 0.0058;
Lr = 0.0058;
Lm = 0.0055;
P  = 2;
J  = 0.02;
B  = 0.001;

sigma = 1 - (Lm^2)/(Ls*Lr);

%% ---------------- CONTROL LIMITS ----------------
Vdc = 300;
Vmax = Vdc/sqrt(3);
Imax = 10;

%% ---------------- REFERENCES ----------------
w_ref = 1500*(2*pi)/(60);          % rad/s
id_ref = 2;           % Flux current

%% ---------------- PI GAINS (TUNED FOR STABILITY) ----------------
Kp_w = 0.5;  Ki_w = 20;
Kp_i = 5;    Ki_i = 300;

%% ---------------- INITIAL STATES ----------------
id = 0; iq = 0;
psi_r = 0.1;
w_m = 0;

int_w  = 0;
int_id = 0;
int_iq = 0;

theta_e = 0;

%% ---------------- FAULT SETUP ----------------
t_fault = 1;
fault_gain = 1;

%% ---------------- DATA STORAGE ----------------
w_log = zeros(1,N);
id_log = zeros(1,N);
iq_log = zeros(1,N);

%% ================= MAIN LOOP =================
for k = 1:N

    %% ---- FAULT INJECTION ----
    if t(k) > t_fault
        fault_gain = 1.5;     % Resistance increase
    end
    Rs_eff = Rs * fault_gain;

    %% ---- SPEED CONTROLLER ----
    w_err = w_ref - w_m;
    int_w = int_w + Ki_w*Ts*w_err;
    int_w = max(min(int_w, Imax), -Imax);

    iq_ref = Kp_w*w_err + int_w;
    iq_ref = max(min(iq_ref, Imax), -Imax);

    %% ---- CURRENT CONTROLLERS ----
    err_id = id_ref - id;
    err_iq = iq_ref - iq;

    int_id = int_id + Ki_i*Ts*err_id;
    int_iq = int_iq + Ki_i*Ts*err_iq;

    vd = Kp_i*err_id + int_id;
    vq = Kp_i*err_iq + int_iq;

    %% ---- VOLTAGE SATURATION ----
    Vmag = sqrt(vd^2 + vq^2);
    if Vmag > Vmax
        vd = vd * Vmax / Vmag;
        vq = vq * Vmax / Vmag;
    end

    %% ---- ELECTRICAL DYNAMICS ----
    did = (vd - Rs_eff*id + w_m*Ls*iq) / Ls;
    diq = (vq - Rs_eff*iq - w_m*Ls*id) / Ls;

    id = id + Ts*did;
    iq = iq + Ts*diq;

    %% ---- TORQUE & MECHANICS ----
    Te = (3/2)*P*(Lm/Lr)*psi_r*iq;
    dw = (Te - B*w_m) / J;
    w_m = w_m + Ts*dw;

    %% ---- FLUX OBSERVER (NORMALIZED) ----
    psi_r = psi_r + Ts*(Lm/Lr*id - psi_r*Rr/Lr);
    psi_r = max(psi_r, 0.05);  % Prevent collapse

    %% ---- ELECTRICAL ANGLE ----
    theta_e = theta_e + Ts*w_m*P;

    %% ---- LOGGING ----
    w_log(k) = w_m;
    id_log(k) = id;
    iq_log(k) = iq;

end

%% ================= PLOTS =================
figure;
subplot(3,1,1)
plot(t, w_log*60/(2*pi), 'LineWidth',1.5)
ylabel('Speed (RPM)')
grid on
subplot(3,1,2)
plot(t, id_log, 'LineWidth',1.5)
ylabel('i_d (A)')
grid on
subplot(3,1,3)
plot(t, iq_log, 'LineWidth',1.5)
ylabel('i_q (A)')
xlabel('Time (s)')
grid on



% figure('Name','Speed');
% plot(t, w_log*60/(2*pi), 'LineWidth',1.5);
% ylabel('Speed (RPM)');
% xlabel('Time (s)');
% grid on; hold on;
% 
% figure('Name','Id');
% plot(t, id_log, 'LineWidth',1.5);
% ylabel('i_d (A)');
% xlabel('Time (s)');
% grid on; hold on;
% 
% figure('Name','Iq');
% plot(t, iq_log, 'LineWidth',1.5);
% ylabel('i_q (A)');
% xlabel('Time (s)');
% grid on; hold on;