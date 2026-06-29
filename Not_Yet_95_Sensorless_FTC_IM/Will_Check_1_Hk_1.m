%% STEP 2 - FOC Integrated Plant: IM in synchronous dq, fed by 2-level VSI (SPWM)
clear; clc; close all;

%% -------------------- Motor Parametreleri (1.5 kW, 4-kutup) [1] ----------
p = 2;              % kutup çifti [3]
Rs = 4.10;          % [ohm] [1]
Rr = 3.90;          % [ohm] [1]
Lls = 8.0e-3;       % [H] [1]
Llr = 8.0e-3;       % [H] [1]
Lm = 0.135;         % [H] [1]
J = 0.008;          % [kg.m^2] [1]
B = 5e-4;           % [N.m.s] [1]
Ls = Lls + Lm;      % [1]
Lr = Llr + Lm;      % [1]

% dq akımlarından dq akılarına geçiş matrisi [1]
Lmat = [ Ls 0 Lm 0; 0 Ls 0 Lm; Lm 0 Lr 0; 0 Lm 0 Lr ];

%% -------------------- Simulasyon Ayarları [4] ---------------------------
Vdc = 560;          % [V] DC bara voltajı [5]
fsw = 10e3;         % Anahtarlama frekansı [5]
Ts = 1e-4;          % Örnekleme zamanı (100 us) [4]
Tend = 1.5;         % Simülasyon süresi
t = (0:Ts:Tend).';
N = numel(t);

% Durum Değişkenleri: x = [ids iqs idr iqr wm]^T [4]
x = zeros(5,1);

% Kayıt Değişkenleri [4, 6]
rpm = zeros(N,1); Te = zeros(N,1);
ids_log = zeros(N,1); iqs_log = zeros(N,1);
va_log = zeros(N,1); vb_log = zeros(N,1); vc_log = zeros(N,1);

%% -------------------- FOC Kontrolcü Parametreleri (Harici Bilgi) ---------
% Bu bölüm kaynaklarda bulunmamaktadır; genel kontrol tasarımıdır.
w_ref = 1000 * (2*pi/60); % 1000 RPM referans hız
ids_ref = 2.5;            % Akı akımı referansı (sabit)

% PI Kazançları
Kp_w = 0.8; Ki_w = 15;    % Hız kontrolcü
Kp_i = 12;  Ki_i = 150;   % Akım kontrolcü

% Entegratör başlangıç değerleri
err_w_int = 0; err_id_int = 0; err_iq_int = 0;
theta_e = 0;              % Elektriksel açı

%% -------------------- Ana Döngü [7] -------------------------------------
Tl = 10; % Yük torku [4]

for k = 1:N
    % --- 1. Geri Besleme Değerlerini Al ---
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_k = x(5);
    wr_e = p * wm_k; % Rotor elektriksel hızı [8]

    % --- 2. Hız Kontrolcü (Dış Döngü - Harici Bilgi) ---
    err_w = w_ref - wm_k;
    err_w_int = err_w_int + err_w * Ts;
    iqs_ref = Kp_w * err_w + Ki_w * err_w_int; % Moment akımı referansı

    % --- 3. Kayma ve Açı Hesaplama (Indirect FOC - Harici Bilgi) ---
    % w_sl = (Rr/Lr) * (iqs_ref / ids_ref)
    wsl = (Rr / Lr) * (iqs_ref / ids_ref); 
    we = wr_e + wsl; % Senkron hız
    theta_e = theta_e + we * Ts; % Akı açısı entegrasyonu
    theta_e = mod(theta_e, 2*pi);

    % --- 4. Akım Kontrolcüleri (İç Döngü - Harici Bilgi) ---
    % D-ekseni (Akı kontrolü)
    err_id = ids_ref - ids_k;
    err_id_int = err_id_int + err_id * Ts;
    vds_ref = Kp_i * err_id + Ki_i * err_id_int;

    % Q-ekseni (Moment kontrolü)
    err_iq = iqs_ref - iqs_k;
    err_iq_int = err_iq_int + err_iq * Ts;
    vqs_ref = Kp_i * err_iq + Ki_i * err_iq_int;

    % --- 5. Ters Dönüşüm (dq -> abc) ve SPWM Modülasyonu ---
    [va_ref, vb_ref, vc_ref] = dq_to_abc(vds_ref, vqs_ref, theta_e);
    
    % Görev oranları (Duty ratios) ve sınırlama [9]
    da = min(max(0.5 + va_ref/Vdc, 0), 1);
    db = min(max(0.5 + vb_ref/Vdc, 0), 1);
    dc = min(max(0.5 + vc_ref/Vdc, 0), 1);

    % Evirici faz voltajları [9]
    va_k = (2*da - 1)*(Vdc/2);
    vb_k = (2*db - 1)*(Vdc/2);
    vc_k = (2*dc - 1)*(Vdc/2);
    v_cm = (va_k + vb_k + vc_k)/3; % Ortak mod voltajı [8]
    va_k = va_k - v_cm; vb_k = vb_k - v_cm; vc_k = vc_k - v_cm;

    % --- 6. Motor Dinamik Modeli (Plant) [2] ---
    [vds, vqs] = abc_to_dq(va_k, vb_k, vc_k, theta_e); % [10]

    % Akılar ve Tork Hesabı [2]
    i4 = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat * i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);
    
    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k); % Elektromanyetik tork [2]

    % Türevler [2]
    dpsi = zeros(4,1);
    dpsi(1) = vds - Rs*ids_k + we*psiqs;
    dpsi(2) = vqs - Rs*iqs_k - we*psids;
    dpsi(3) = - Rr*idr_k + wsl*psiqr;
    dpsi(4) = - Rr*iqr_k - wsl*psidr;

    di = Lmat \ dpsi;
    dwm = (Te_k - Tl - B*wm_k)/J;

    % İntegrasyon [2]
    x(1:4) = x(1:4) + Ts*di;
    x(5) = x(5) + Ts*dwm;

    % Kayıt [6]
    rpm(k) = wm_k*60/(2*pi);
    Te(k) = Te_k;
    ids_log(k) = ids_k; iqs_log(k) = iqs_k;
    va_log(k) = va_k; vb_log(k) = vb_k; vc_log(k) = vc_k;
end

%% -------------------- Grafikler [6, 11] ---------------------------------
figure; subplot(2,1,1);
plot(t, rpm, 'b', 'LineWidth', 1.5); grid on; ylabel('Hız [RPM]');
title('FOC Kontrollü Asenkron Motor');
subplot(2,1,2);
plot(t, Te, 'r', 'LineWidth', 1.5); grid on; ylabel('Tork [Nm]');
xlabel('Zaman [s]');

figure; plot(t, ids_log, t, iqs_log, 'LineWidth', 1.2); grid on;
legend('i_{ds}','i_{qs}'); title('dq Eksen Akımları');

%% -------------------- Yerel Fonksiyonlar [10] ---------------------------
function [vd, vq] = abc_to_dq(va, vb, vc, theta)
    alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    beta = (2/3)*((sqrt(3)/2)*(vb - vc));
    vd = alpha*cos(theta) + beta*sin(theta);
    vq = -alpha*sin(theta) + beta*cos(theta);
end

function [va, vb, vc] = dq_to_abc(vd, vq, theta)
    % Ters Park ve Ters Clarke (Harici Bilgi)
    alpha = vd*cos(theta) - vq*sin(theta);
    beta  = vd*sin(theta) + vq*cos(theta);
    va = alpha;
    vb = -0.5*alpha + (sqrt(3)/2)*beta;
    vc = -0.5*alpha - (sqrt(3)/2)*beta;
end