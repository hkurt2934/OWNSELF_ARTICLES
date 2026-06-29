clc; clear; close all;
% Devam edilmesi gereken: Chat GPT Hkurt2934 deki Important önerisini buraya
% uygulanması.


% Makale gücünü belirleyen asıl faktör (kritik nokta)
% FOC tipi tek başına makaleyi güçlü yapmaz.
% Aşağıdaki soruların “evet” cevabı varsa makale güçlüdür:
% Arıza altında akı yönlendirme bozulmadan devam ediyor mu?
% Observer arıza koşulunda kararlılığını koruyor mu?
% DC-link sınırlaması + SVPWM birlikte ele alınmış mı?
% Anti-windup ve saturasyon kontrol mimarisine entegre mi?
% Torque ripple ve hız salınımı nicel olarak gösterilmiş mi?
% Bu sorular Direct FOC ile çok daha iyi cevaplanır.

% Sadece IFOC	❌ Zayıf
% IFOC + küçük iyileştirme	⚠️ Orta
% DFOC + observer	✅ Güçlü
% DFOC + observer + FTC	🔥 Çok güçlü
% DFOC + FTC + DC-link + SVPWM	🚀 Q1 adayı

% Güçlü, uzun ömürlü ve Q1 hedefli bir makale için:
% Direct FOC tabanlı, observer destekli ve arıza toleranslı bir yapı seçmelisin.

% Neden Direkt SVPWM yok:
% SVPWM is employed for inverter control.
% To focus on the machine dynamics and sensorless control performance, an average-value inverter model is used in simulations.”
% Bu ifade:
% Literatürde standarttır
% IEEE makalelerinde yüzlerce örneği vardır
% Hakem için tamamen yeterlidir

% Ortalama inverter kullanmak neyi BOZMAZ?
% Bozmaz:
% Rotor akı yönelimi
% Id–Iq ayrışması
% MRAS hız kestirimi
% Speed loop stabilitesi
% Load torque response
% Low-speed instability analizi

% Bozmaz:
% “Sensorless” iddiasını
% “Direct FOC” tanımını
% Akı gözleyici geçerliliğini

% Electrical machines hakemi neye bakar?
% Hakem şunları görmek ister:
% ✔ 1. Akı düzgün mü yönlenmiş?
% ψrq≈0
% ψrd sabit mi?
% ✔ 2. Id–Iq ayrışması var mı?
% Yük değişiminde Id sabit mi?
% Iq torka mı gidiyor?
% ✔ 3. MRAS stabil mi?
% Düşük hızda sapma var mı?
% Yük adımında overshoot makul mü?
% ✔ 4. Speed loop düzgün mü?
% Overshoot kabul edilebilir mi?
% Settling time mantıklı mı?
% Bunların hiçbiri anahtarlamalı SVPWM gerektirmez.

% Şu tür cümleleri IEEE makalelerinde sürekli görürsün:
% “An average inverter model is used to reduce computational burden.”
% “Switching dynamics are neglected to emphasize control performance.”
% “The proposed sensorless scheme is validated using an average-value VSI.”
% Bu ifadeler red sebebi değildir.

% İyi bir electrical machines makalesi için:
% ✔ Ortalama inverter
% ✔ SVPWM matematiksel sınırları
% ✔ Temiz FOC
% ✔ Stabil MRAS
% ✔ Doğru blok diyagram
% Bu yapı akademik olarak tamamen doğrudur ve güvenlidir.

%% ===============================
% 1. Motor Parameters
%% ===============================
Rs = 1.405;      % Stator resistance (ohm)
Rr = 1.395;      % Rotor resistance (ohm)
Ls = 0.0058;     % Stator inductance (H)
Lr = 0.0058;     % Rotor inductance (H)
Lm = 0.0055;     % Magnetizing inductance (H)
pp  = 2;          % Pole pairs
J  = 0.02;       % Inertia (kg.m^2)
B  = 0.001;      % Friction coefficient

Tr = Lr / Rr;    % Rotor time constant

%% ===============================
% 2. Simulation Parameters
%% ===============================
Ts   = 1e-4;
Tsim = 3;
N    = Tsim/Ts;
t    = (0:N-1)*Ts;

Vdc = 560;
Vmax = Vdc/sqrt(3);
% iq_max = 1.2 * In;   % 120% rated current
TL  = 10;

%% ===============================
% 3. Control Parameters
%% ===============================
Kpw = 1.5;    Kiw = 50;
Kpd = 50;   Kid = 500;
Kpq = 50;   Kiq = 500;
Kaw_id = 50/Kpd;    % or 1/Kpq for q-axis
Kaw_iq = 1/Kpq;

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
w_e = 0;

%% ===============================
% 5. Logging
%% ===============================
wr_log = zeros(1,N); wr_hat_log = zeros(1,N);
ids_log=zeros(N,1); iqs_log=zeros(N,1);
Te_log=zeros(N,1);

%% ===============================
% References
%% ===============================
psi_r_ref = 0.8;    % Flux Reference to calculate id reference
w_ref = 100;        % rad/s
ids_ref = psi_r_ref/Lm;

psi_r = 0;
psi_r_log = zeros(N,1);

%% ===============================
% 6. Main Loop
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

    % SVPWM voltage magnitude
    Vref = sqrt(vds_ref^2 + vqs_ref^2);

    % DC-link constraint (SVPWM)
    if Vref > Vmax
        vds = vds_ref * Vmax / Vref;
        vqs = vqs_ref * Vmax / Vref;
    else
        vds = vds_ref;
        vqs = vqs_ref;
    end

    % Anti-windup (back-calculation)
    int_id = int_id + (vds - vds_ref)*Kaw_id*Ts;
    int_iq = int_iq + (vqs - vqs_ref)*Kaw_iq*Ts;

    %% -------- Inverse Park (dq → αβ) ---
    v_alpha =  vds*cos(theta_e) - vqs*sin(theta_e);
    v_beta  =  vds*sin(theta_e) + vqs*cos(theta_e);

    %% -------- Inverse Clarke (αβ → abc)
    va = v_alpha;
    vb = -0.5*v_alpha + sqrt(3)/2*v_beta;
    vc = -0.5*v_alpha - sqrt(3)/2*v_beta;

    %% ==================================
    %  Induction Motor dq Model
    %% ==================================
    dpsi_r = (Lm/Tr)*ids - (1/Tr)*psi_r;
    psi_r  = psi_r + Ts*dpsi_r;
    
    psi_r_eps = 1e-4;
    psi_r_eff = max(abs(psi_r), psi_r_eps) * sign(psi_r + psi_r_eps);
    
    wsl = (Lm/Tr) * (iqs / psi_r_eff);     % electrical rad/s

    we = wr + wsl;                         % if wr is rotor electrical speed
    % If wr is mechanical speed (rad/s), use: we = p*wr + wsl;

    % --- Electromagnetic torque (your requested formula) ---
    Te = (3/2)*pp*(Lm/Lr)*psi_r*iqs;
    
    dwr = (Te - TL - B*wr)/J;
    wr  = wr + dwr*Ts;



    % dpsi_rd = (Lm/Tr)*ids - psi_rd/Tr + wr*psi_rq;
    % dpsi_rq = (Lm/Tr)*iqs - psi_rq/Tr - wr*psi_rd;
    % 
    % psi_rd = psi_rd + Ts*dpsi_rd;
    % psi_rq = psi_rq + Ts*dpsi_rq;

    dids = (vds - Rs*ids + wr*Ls*iqs)/Ls;
    diqs = (vqs - Rs*iqs - wr*Ls*ids)/Ls;

    ids = ids + Ts*dids;
    iqs = iqs + Ts*diqs;

    %% -------- Inverse Park (dq → αβ currents)
    i_alpha =  ids*cos(theta_e) - iqs*sin(theta_e);
    i_beta  =  ids*sin(theta_e) + iqs*cos(theta_e);

    %% -------- Inverse Clarke (αβ → abc currents)
    ia = i_alpha;
    ib = -0.5*i_alpha + sqrt(3)/2*i_beta;
    ic = -0.5*i_alpha - sqrt(3)/2*i_beta;

    %% ==================================
    %  MRAS Speed Estimator
    %% ==================================
    % Clarke voltages
    v_alpha_m = va;
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
    dpsi_ad_alpha = -psi_ad(1)/Tr + w_e*psi_ad(2) + (Lm/Tr)*i_alpha_m;
    dpsi_ad_beta  = -psi_ad(2)/Tr - w_e*psi_ad(1) + (Lm/Tr)*i_beta_m;

    psi_ad(1) = psi_ad(1) + Ts*dpsi_ad_alpha;
    psi_ad(2) = psi_ad(2) + Ts*dpsi_ad_beta;

    % MRAS error
    e_mras = psi_ref(1)*psi_ad(2) - psi_ref(2)*psi_ad(1);

    % Speed adaptation
    int_mras = int_mras + e_mras*Ts;
    w_e = Kp_mras*e_mras + Ki_mras*int_mras;

    % Electrical angle
    theta_e = theta_e + w_e*Ts;

    %% -------- Logging
    wr_log(k) = wr;
    wr_hat_log(k) = w_e;
    ids_log(k)=ids;
    iqs_log(k)=iqs;
    psi_r_log(k) = psi_r;
    Te_log(k)=Te;

end

%% ===============================
% 7. Results
%% ===============================

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

