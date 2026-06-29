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

Vdc = 560;
TL  = 10;

%% ===============================
% 3. Control Parameters
%% ===============================
Kpw = 2;    Kiw = 50;
Kpd = 30;   Kid = 500;
Kpq = 30;   Kiq = 500;

Kp_mras = 200;
Ki_mras = 5000;

%% ===============================
% 4. Initialization
%% ===============================
ids = 0; iqs = 0;
psi_rd = 0; psi_rq = 0;
wr = 0;

theta_e = 0;

% Controllers
int_w = 0; int_id = 0; int_iq = 0; int_mras = 0;

% MRAS
psi_ref = [0;0];
psi_ad  = [0;0];
omega_hat = 0;
theta_hat = 0;

%% ===============================
% 5. Logging
%% ===============================
wr_log = zeros(1,N);
wr_hat_log = zeros(1,N);

%% ===============================
% References
%% ===============================
psi_r_ref = 0.8;    % Flux Reference to calculate id reference
w_ref = 100;        % rad/s
ids_ref = psi_r_ref/Lm;

%% ===============================
% 6. Main Loop
%% ===============================
for k = 1:N

    %% -------- Speed Controller --------
    omega_hat = wr;
    ew = w_ref - omega_hat;
    int_w = int_w + ew*Ts;
    iqs_ref = Kpw*ew + Kiw*int_w;

    %% -------- Current Controllers -----
    ed = ids_ref - ids;
    int_id = int_id + ed*Ts;
    vds = Kpd*ed + Kid*int_id;

    eq = iqs_ref - iqs;
    int_iq = int_iq + eq*Ts;
    vqs = Kpq*eq + Kiq*int_iq;

    %% -------- Inverse Park (dq → αβ) ---
    v_alpha =  vds*cos(theta_hat) - vqs*sin(theta_hat);
    v_beta  =  vds*sin(theta_hat) + vqs*cos(theta_hat);

    %% -------- Inverse Clarke (αβ → abc)
    va = v_alpha;
    vb = -0.5*v_alpha + sqrt(3)/2*v_beta;
    vc = -0.5*v_alpha - sqrt(3)/2*v_beta;

    %% ==================================
    %  Induction Motor dq Model
    %% ==================================
    Te = (3/2)*P*(Lm/Lr)*(psi_rd*iqs - psi_rq*ids); % ??? Neden iqs ve ids direkt alındı burada. öncesinde hesaplama yok.

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

    %% -------- Inverse Park (dq → αβ currents)
    i_alpha =  ids*cos(theta_hat) - iqs*sin(theta_hat);
    i_beta  =  ids*sin(theta_hat) + iqs*cos(theta_hat);

    %% -------- Inverse Clarke (αβ → abc currents)
    ia = i_alpha;
    ib = -0.5*i_alpha + sqrt(3)/2*i_beta;
    ic = -0.5*i_alpha - sqrt(3)/2*i_beta;

    %% ==================================
    %  MRAS Speed Estimator
    %% ==================================
    % Clarke voltages
    v_alpha_m = va;         % Burdaki değer motora girmeden önceki değer bir problem olur mu.
    v_beta_m  = (va + 2*vb)/sqrt(3);

    % Clarke currents
    i_alpha_m = ia;
    i_beta_m  = (ia + 2*ib)/sqrt(3);

    % Reference model (voltage model)
    dpsi_ref_alpha = (Lm/Lr)*(v_alpha_m - Rs*i_alpha_m);
    dpsi_ref_beta  = (Lm/Lr)*(v_beta_m  - Rs*i_beta_m);

    psi_ref(1) = psi_ref(1) + Ts*dpsi_ref_alpha;
    psi_ref(2) = psi_ref(2) + Ts*dpsi_ref_beta;

    % Adaptive model (current model)
    dpsi_ad_alpha = -psi_ad(1)/Tr + omega_hat*psi_ad(2) + (Lm/Tr)*i_alpha_m;
    dpsi_ad_beta  = -psi_ad(2)/Tr - omega_hat*psi_ad(1) + (Lm/Tr)*i_beta_m;

    psi_ad(1) = psi_ad(1) + Ts*dpsi_ad_alpha;
    psi_ad(2) = psi_ad(2) + Ts*dpsi_ad_beta;

    % MRAS error
    eps = psi_ref(1)*psi_ad(2) - psi_ref(2)*psi_ad(1);

    % Speed adaptation
    int_mras = int_mras + eps*Ts;
    omega_hat = Kp_mras*eps + Ki_mras*int_mras;

    % Electrical angle
    theta_hat = theta_hat + omega_hat*Ts;

    %% -------- Logging
    wr_log(k) = wr;
    wr_hat_log(k) = omega_hat;

end

%% ===============================
% 7. Results
%% ===============================
figure;
plot(t,wr_log,'b','LineWidth',1.5); hold on;
plot(t,wr_hat_log,'r--','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Speed (rad/s)');
legend('Actual','Estimated');
title('Sensorless IM FOC with MRAS');
