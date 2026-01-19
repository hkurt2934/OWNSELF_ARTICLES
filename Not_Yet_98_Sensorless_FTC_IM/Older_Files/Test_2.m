%% COMPREHENSIVE COMPARISON: FTC vs NON-FTC
clear; clc; close all;

%% ---------------- SIMULATION PARAMETERS ----------------
Ts = 1e-4; Tsim = 3; 
N = round(Tsim/Ts);
t = (0:N-1)*Ts;

%% ---------------- MOTOR PARAMETERS ----------------
Rs = 1.405; Rr = 1.395;
Ls = 0.0058; Lr = 0.0058;
Lm = 0.0055; P = 2;
J = 0.02; B = 0.001;

%% ---------------- LIMITS ----------------
Vdc = 300; Vmax = Vdc/sqrt(3);
Imax = 10;

%% ---------------- REFERENCES ----------------
w_ref = 1500*2*pi/60;
id_ref = 2;

%% ---------------- PI GAINS ----------------
Kp_w = 0.5; Ki_w = 20;
Kp_i = 5; Ki_i = 300;

%% ---------------- FAULT ----------------
t_fault = 1;
Rs_fault = 1.5 * Rs;

%% ---------------- CASE DEFINITIONS ----------------
% [fault, FTC, sensorless]
cases = [
    0 0 0;   % C1: Healthy Sensor-based
    0 0 1;   % C2: Healthy Sensorless
    1 0 1;   % C3: Faulty Sensorless (NO FTC)
    1 1 1    % C4: Faulty Sensorless (FTC)
];

labels = {
    'Healthy Sensor-Based';
    'Healthy Sensorless';
    'Faulty Sensorless (No FTC)';
    'Faulty Sensorless (FTC)'
};

%% ---------------- STORAGE ----------------
w_all = zeros(4,N);
iq_all = zeros(4,N);

%% ================= MAIN LOOP =================
for c = 1:4

    fault = cases(c,1);
    FTC   = cases(c,2);
    sens  = cases(c,3);

    id=0; iq=0; psi_r=0.1; w_m=0; w_est=0;
    int_w=0; int_id=0; int_iq=0;

    for k = 1:N

        %% PLANT RESISTANCE
        if fault && t(k)>=t_fault
            Rs_real = Rs_fault;
        else
            Rs_real = Rs;
        end

        %% CONTROLLER RESISTANCE (FTC EFFECT)
        if FTC
            Rs_ctrl = Rs_real;
        else
            Rs_ctrl = Rs;
        end

        %% SPEED FEEDBACK
        if sens
            w_fb = w_est;
        else
            w_fb = w_m;
        end

        %% SPEED PI
        w_err = w_ref - w_fb;
        int_w = max(min(int_w + Ki_w*Ts*w_err, Imax), -Imax);
        iq_ref = max(min(Kp_w*w_err + int_w, Imax), -Imax);

        %% CURRENT PI
        err_id = id_ref - id;
        err_iq = iq_ref - iq;
        int_id = int_id + Ki_i*Ts*err_id;
        int_iq = int_iq + Ki_i*Ts*err_iq;

        vd = Kp_i*err_id + int_id;
        vq = Kp_i*err_iq + int_iq;

        %% SATURATION
        Vmag = hypot(vd,vq);
        if Vmag>Vmax
            vd=vd*Vmax/Vmag; vq=vq*Vmax/Vmag;
        end

        %% MOTOR MODEL
        did = (vd - Rs_real*id + w_m*Ls*iq)/Ls;
        diq = (vq - Rs_real*iq - w_m*Ls*id)/Ls;
        id = id + Ts*did; iq = iq + Ts*diq;

        Te = (3/2)*P*(Lm/Lr)*psi_r*iq;
        w_m = w_m + Ts*(Te - B*w_m)/J;

        psi_r = max(psi_r + Ts*(Lm/Lr*id - psi_r*Rr/Lr), 0.05);
        w_est = w_est + Ts*(Te - B*w_est)/J;

        %% LOGGING
        w_all(c,k) = w_m;
        iq_all(c,k) = iq;
    end
end

%% ================= PLOTS =================
figure('Name','All');
plot(t, w_all'*60/(2*pi), 'LineWidth',1.5);
legend(labels,'Location','best');
ylabel('Speed (rpm)');
xlabel('Time (s)');
grid on; hold on;

figure('Name','Healthy Sensor-based');
plot(t, w_all'*60/(2*pi), 'LineWidth',1.5);
legend(labels,'Location','best');
ylabel('Speed (rpm)');
xlabel('Time (s)');
grid on; hold on;

figure('Name','All');
plot(t, w_all'*60/(2*pi), 'LineWidth',1.5);
legend(labels,'Location','best');
ylabel('Speed (rpm)');
xlabel('Time (s)');
grid on; hold on;

figure('Name','All');
plot(t, w_all'*60/(2*pi), 'LineWidth',1.5);
legend(labels,'Location','best');
ylabel('Speed (rpm)');
xlabel('Time (s)');
grid on; hold on;
