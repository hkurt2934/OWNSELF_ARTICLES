clc; clear; close all;

%% ===============================
% 1. Motor Parameters
%% ===============================
Rs = 1.405;
Rr = 1.395;
Ls = 0.0058;
Lr = 0.0058;
Lm = 0.0055;
P  = 2;

J = 0.02;
B = 0.001;

Tr = Lr/Rr;

%% ===============================
% 2. Simulation Parameters
%% ===============================
Ts   = 1e-4;
Tsim = 3;
N    = Tsim/Ts;
t    = (0:N-1)*Ts;

TL = 10;

%% ===============================
% 3. Control Parameters
%% ===============================
Kpw = 2;    Kiw = 50;
Kpd = 30;   Kid = 500;
Kpq = 30;   Kiq = 500;

psi_r_ref = 0.8;
w_ref = 100;

%% ===============================
% 4. Initialization
%% ===============================
ids = 0; iqs = 0;
psi_rd = 0; psi_rq = 0;
wr = 0;

theta_e = 0;

int_w = 0; int_id = 0; int_iq = 0;

%% ===============================
% 5. Logging
%% ===============================
wr_log = zeros(1,N);

%% ===============================
% 6. Main Loop
%% ===============================
for k = 1:N

    %% -------- Speed Measurement --------
    wr_meas = wr;   % Encoder / speed sensor

    %% -------- Speed Controller --------
    ew = w_ref - wr_meas;
    int_w = int_w + ew*Ts;
    iqs_ref = Kpw*ew + Kiw*int_w;
    ids_ref = psi_r_ref/Lm;

    %% -------- Current Controllers -----
    ed = ids_ref - ids;
    int_id = int_id + ed*Ts;
    vds = Kpd*ed + Kid*int_id;

    eq = iqs_ref - iqs;
    int_iq = int_iq + eq*Ts;
    vqs = Kpq*eq + Kiq*int_iq;

    %% -------- Electrical Angle --------
    theta_e = theta_e + wr_meas*Ts;

    %% -------- Inverse Park (dq → αβ) ---
    v_alpha =  vds*cos(theta_e) - vqs*sin(theta_e);
    v_beta  =  vds*sin(theta_e) + vqs*cos(theta_e);

    %% -------- Inverse Clarke (αβ → abc)
    va = v_alpha;
    vb = -0.5*v_alpha + sqrt(3)/2*v_beta;
    vc = -0.5*v_alpha - sqrt(3)/2*v_beta;

    %% ===============================
    %  Induction Motor dq Model
    %% ===============================
    Te = (3/2)*P*(Lm/Lr)*(psi_rd*iqs - psi_rq*ids);

    dwr = (Te - TL - B*wr)/J;
    wr  = wr + dwr*Ts;

    dpsi_rd = (Lm/Tr)*ids - psi_rd/Tr + wr*psi_rq;
    dpsi_rq = (Lm/Tr)*iqs - psi_rq/Tr - wr*psi_rd;

    psi_rd = psi_rd + Ts*dpsi_rd;
    psi_rq = psi_rq + Ts*dpsi_rq;

    dids = (vds - Rs*ids + wr*Ls*iqs)/Ls;
    diqs = (vqs - Rs*iqs - wr*Ls*ids)/Ls;

    ids = ids + Ts*dids;
    iqs = iqs + Ts*diqs;

    %% -------- Logging --------
    wr_log(k) = wr;

end

%% ===============================
% 7. Plot
%% ===============================
figure;
plot(t,wr_log,'b','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Speed (rad/s)');
title('Sensor-Based FOC of Induction Motor');
