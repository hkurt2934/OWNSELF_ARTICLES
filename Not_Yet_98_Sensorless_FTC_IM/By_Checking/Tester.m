clc; clear;

%% =========================================================
% 1. MOTOR PARAMETRELERİ (IM)
%% =========================================================
Rs = 1.405;
Rr = 1.395;
Ls = 0.0058;
Lr = 0.0058;
Lm = 0.0055;

P  = 2;          % Pole pairs
J  = 0.02;
B  = 0.001;

Tr = Lr/Rr;
sigma = 1 - (Lm^2)/(Ls*Lr);

%% =========================================================
% 2. SİMÜLASYON PARAMETRELERİ
%% =========================================================
Ts   = 1e-4;
Tsim = 3;
N    = Tsim/Ts;
t    = (0:N-1)*Ts;

Vdc  = 560;
Vmax = Vdc/sqrt(3);

TL = 10;

%% =========================================================
% 3. KONTROL PARAMETRELERİ
%% =========================================================
% Speed PI (MRAS speed)
Kpw = 1.5;   Kiw = 50;

% Current PI
Kpd = 40;  Kid = 800;
Kpq = 40;  Kiq = 800;

% MRAS
Kp_mras = 200;
Ki_mras = 5000;

%% =========================================================
% 4. BAŞLANGIÇ DEĞERLERİ
%% =========================================================
ids = 0; iqs = 0;
wr  = 0;

psi_rd = 0; psi_rq = 0;   % only for motor model

psi_s_alpha = 0;
psi_s_beta  = 0;

theta_e = 0;

% Integrators
int_w = 0; int_id = 0; int_iq = 0; int_mras = 0;

% MRAS fluxes
psi_ref = [0;0];
psi_ad  = [0;0];

w_hat = 0;

%% =========================================================
% 5. REFERANSLAR
%% =========================================================
psi_s_ref = 0.8;
ids_ref   = psi_s_ref/Lm;
w_ref     = 100;   % rad/s

%% =========================================================
% 6. LOG
%% =========================================================
wr_log = zeros(N,1);
wr_hat_log = zeros(N,1);
Te_log = zeros(N,1);
ids_log=zeros(N,1); iqs_log=zeros(N,1);

%% =========================================================
% 7. ANA DÖNGÜ
%% =========================================================
for k = 1:N

    %% ---------------------------------
    %  SPEED CONTROLLER (sensorless)
    %% ---------------------------------
    ew = w_ref - w_hat;
    int_w = int_w + ew*Ts;
    iqs_ref = Kpw*ew + Kiw*int_w;

    %% ---------------------------------
    %  CURRENT CONTROLLERS (dq)
    %% ---------------------------------
    ed = ids_ref - ids;
    int_id = int_id + ed*Ts;
    vds_ref = Kpd*ed + Kid*int_id;

    eq = iqs_ref - iqs;
    int_iq = int_iq + eq*Ts;
    vqs_ref = Kpq*eq + Kiq*int_iq;

    % SVPWM saturation
    Vref = sqrt(vds_ref^2 + vqs_ref^2);
    if Vref > Vmax
        vds = vds_ref * Vmax/Vref;
        vqs = vqs_ref * Vmax/Vref;
    else
        vds = vds_ref;
        vqs = vqs_ref;
    end

    %% ---------------------------------
    %  Inverse Park
    %% ---------------------------------
    v_alpha =  vds*cos(theta_e) - vqs*sin(theta_e);
    v_beta  =  vds*sin(theta_e) + vqs*cos(theta_e);

    %% ---------------------------------
    %  Inverse Clarke
    %% ---------------------------------
    va = v_alpha;
    vb = -0.5*v_alpha + sqrt(3)/2*v_beta;
    vc = -0.5*v_alpha - sqrt(3)/2*v_beta;

    %% =====================================================
    %  INDUCTION MOTOR MODEL (SIMULATION ONLY)
    %% =====================================================
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

    %% ---------------------------------
    %  Currents abc → αβ
    %% ---------------------------------
    i_alpha = ids*cos(theta_e) - iqs*sin(theta_e);
    i_beta  = ids*sin(theta_e) + iqs*cos(theta_e);

    ia = i_alpha;
    ib = -0.5*i_alpha + sqrt(3)/2*i_beta;
    ic = -0.5*i_alpha - sqrt(3)/2*i_beta;

    %% =====================================================
    %  STATOR FLUX OBSERVER (DIRECT FOC)
    %% =====================================================
    dpsi_s_alpha = v_alpha - Rs*i_alpha;
    dpsi_s_beta  = v_beta  - Rs*i_beta;

    psi_s_alpha = psi_s_alpha + Ts*dpsi_s_alpha;
    psi_s_beta  = psi_s_beta  + Ts*dpsi_s_beta;

    theta_e = atan2(psi_s_beta, psi_s_alpha);

    %% =====================================================
    %  MRAS SPEED ESTIMATOR (NO ANGLE)
    %% =====================================================
    dpsi_ref_alpha = (Lm/Lr)*(v_alpha - Rs*i_alpha);
    dpsi_ref_beta  = (Lm/Lr)*(v_beta  - Rs*i_beta);

    psi_ref(1) = psi_ref(1) + Ts*dpsi_ref_alpha;
    psi_ref(2) = psi_ref(2) + Ts*dpsi_ref_beta;

    dpsi_ad_alpha = -psi_ad(1)/Tr + w_hat*psi_ad(2) + (Lm/Tr)*i_alpha;
    dpsi_ad_beta  = -psi_ad(2)/Tr - w_hat*psi_ad(1) + (Lm/Tr)*i_beta;

    psi_ad(1) = psi_ad(1) + Ts*dpsi_ad_alpha;
    psi_ad(2) = psi_ad(2) + Ts*dpsi_ad_beta;

    e_mras = psi_ref(1)*psi_ad(2) - psi_ref(2)*psi_ad(1);
    int_mras = int_mras + e_mras*Ts;

    w_hat = Kp_mras*e_mras + Ki_mras*int_mras;

    %% ---------------------------------
    %  LOG
    %% ---------------------------------
    wr_log(k)     = wr;
    w_hat_log(k) = w_hat;
    Te_log(k)     = Te;
    ids_log(k)=ids;
    iqs_log(k)=iqs;

end
wr_rpm = wr_log*60/(2*pi);
wr_ref_rpm = w_ref*60/(2*pi);
figure('Name','Speed');
subplot(2,2,1);
plot(t, w_ref*ones(size(t)),'LineWidth',1.2); hold on;
plot(t, wr_rpm,'LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('Speed [rpm]');
legend('Reference','Actual'); title('FOC Speed Control (Healthy, SPWM avg)');

% figure('Name','Torque');
subplot(2,2,2);
plot(t, Te_log,'LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('T_e [N*m]'); title('Electromagnetic torque');

% figure('Name','dq currents');
subplot(2,2,3);
plot(t, ids_log,'LineWidth',1.1); hold on;
plot(t, iqs_log,'LineWidth',1.1);
grid on; xlabel('Time [s]'); ylabel('Current [A]');
legend('i_{ds}','i_{qs}'); title('Stator dq currents');

% figure('Name','First');
subplot(2,2,4);
plot(t,wr_log,'b','LineWidth',1.5); hold on;
plot(t,wr_hat_log,'r--','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Speed (rad/s)');
legend('Actual','Estimated');
title('Sensorless IM FOC with MRAS');