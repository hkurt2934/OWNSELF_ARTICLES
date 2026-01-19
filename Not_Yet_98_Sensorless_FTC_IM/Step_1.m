%% 1.5 kW Induction Motor fed by 50 Hz 3-phase Voltage Source (Ideal Sinusoidal)
% - Synchronous dq IM model
% - Balanced 3-phase sinusoidal source at fixed 50 Hz
% - Constant load torque (choose rated ~9.55 N.m or heavy 50 N.m)
% - Outputs: speed, torque, dq currents

clear; clc; close all;

%% ----------------------- Target / Ratings --------------------------------
P_rated = 1.5e3;              % [W]
n_sync  = 1500;               % [rpm] synchronous at 50 Hz, 4-pole
w_sync  = n_sync*2*pi/60;     % [rad/s]
T_rated = P_rated / w_sync;   % ~9.55 N.m

%% ------------------------------ Load torque ------------------------------
% Choose ONE:
Tl = T_rated;              % rated-load run (~9.55 N.m)
% Tl = 50;                 % heavy-load stress test (will likely stall for 1.5 kW)

%% ------------------------------ Time -------------------------------
f   = 50;           % [Hz]
Ts   = 1e-4;               % [s]
Tend = 3;                % [s]
t = (0:Ts:Tend).';
N = numel(t);

%% ----------------------- Motor Parameters (typical 1.5 kW set) -----------
% 3-phase squirrel-cage IM, ~1.5 kW class, 4-pole (2 pole-pairs)
p   = 2;            % pole pairs
Rs  = 4.10;         % stator resistance [ohm]
Rr  = 3.90;         % rotor resistance [ohm]
Lls = 8.0e-3;       % stator leakage inductance [H]
Llr = 8.0e-3;       % rotor leakage inductance [H]
Lm  = 0.135;        % magnetizing inductance [H]
J   = 0.008;        % inertia [kg.m^2] (small motor)
B   = 5e-4;         % viscous friction [N.m.s]

Ls  = Lls + Lm;
Lr  = Llr + Lm;

Lsigma = Ls - (Lm^2)/Lr; % leakage inductance (sigma)
Tr = Lr/Rr;              % rotor time constant
theta_e = 0;             % electrical angle

% Inductance matrix for dq currents -> dq fluxes (power-invariant form)
% psi = Lmat * i, where i = [ids iqs idr iqr]^T
Lmat = [ Ls   0   Lm   0;
          0  Ls    0  Lm;
         Lm   0   Lr   0;
          0  Lm    0  Lr ];

%% -------------------------- Supply --------------------------
we  = 2*pi*f;              % [rad/s] electrical angular frequency
Vll_rms = 400;             % [V] line-line RMS (common grid/inverter bus)
Vph_rms = Vll_rms/sqrt(3);
Vph_pk  = sqrt(2)*Vph_rms;

%% -------------------- Simulation ----------------------------------------
% States: x = [ids iqs idr iqr wm]^T
x = zeros(5,1);

% Logs
ids = zeros(N,1); iqs = zeros(N,1);
idr = zeros(N,1); iqr = zeros(N,1);
wm  = zeros(N,1);
Te  = zeros(N,1);

vds_log = zeros(N,1); vqs_log = zeros(N,1);

%% -------------------- Main loop -----------------------------------------
for k = 1:N
    tk = t(k);

    % Synchronous electrical angle for the dq reference frame
    theta_e = we * tk;

    % Balanced 3-phase stator voltages
    va = Vph_pk*sin(theta_e);
    vb = Vph_pk*sin(theta_e - 2*pi/3);
    vc = Vph_pk*sin(theta_e + 2*pi/3);

    % abc -> dq (power-invariant)
    [vds, vqs] = abc_to_dq(va, vb, vc, theta_e);

    % Unpack states
    ids_k = x(1); iqs_k = x(2);
    idr_k = x(3); iqr_k = x(4);
    wm_k  = x(5);

    % Electrical rotor speed and slip frequency
    wr_e = p * wm_k;          % electrical rotor speed [rad/s]
    wsl  = we - wr_e;         % slip electrical speed [rad/s]

    % Fluxes from currents
    i4  = [ids_k; iqs_k; idr_k; iqr_k];
    psi = Lmat * i4;
    psids = psi(1); psiqs = psi(2);
    psidr = psi(3); psiqr = psi(4);

    % Torque
    Te_k = (3/2)*p*(psids*iqs_k - psiqs*ids_k);

    % dq electrical dynamics (in synchronous frame)
    % Stator:
    % vds = Rs*ids + d(psids)/dt - we*psiqs
    % vqs = Rs*iqs + d(psiqs)/dt + we*psids
    % Rotor (shorted):
    % 0 = Rr*idr + d(psidr)/dt - wsl*psiqr
    % 0 = Rr*iqr + d(psiqr)/dt + wsl*psidr
    %
    % Rearranged for d(psi)/dt:
    dpsi = zeros(4,1);
    dpsi(1) = vds - Rs*ids_k + we*psiqs;
    dpsi(2) = vqs - Rs*iqs_k - we*psids;
    dpsi(3) =      - Rr*idr_k + wsl*psiqr;
    dpsi(4) =      - Rr*iqr_k - wsl*psidr;

    % Convert dpsi -> di via Lmat * di = dpsi
    di = Lmat \ dpsi;

    % Mechanical
    dwm = (Te_k - Tl - B*wm_k)/J;

    % Integrate (Euler)
    x(1:4) = x(1:4) + Ts*di;
    x(5)   = x(5)   + Ts*dwm;

    % Log
    ids(k)=ids_k; iqs(k)=iqs_k; idr(k)=idr_k; iqr(k)=iqr_k;
    wm(k)=wm_k; Te(k)=Te_k;
    vds_log(k)=vds; vqs_log(k)=vqs;
end
%% -------------------- Plots ---------------------------------------------
rpm = wm*60/(2*pi);

figure; plot(t, rpm, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Speed [rpm]');
title(sprintf('Healthy IM (Open-loop 50 Hz Voltage), Load = %.2f N·m', Tl));

figure; plot(t, Te, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Torque [N·m]');
title('Electromagnetic Torque');

figure; plot(t, ids, 'LineWidth', 1.2); hold on;
plot(t, iqs, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Current [A]');
title('Stator Currents in synchronous dq');
legend('i_{ds}','i_{qs}');

%% -------------------- Local function ------------------------------------
function [vd, vq] = abc_to_dq(va, vb, vc, theta)
    % Power-invariant Clarke
    alpha = (2/3)*(va - 0.5*vb - 0.5*vc);
    beta  = (2/3)*((sqrt(3)/2)*(vb - vc));

    % Park
    c = cos(theta); s = sin(theta);
    vd =  alpha*c + beta*s;
    vq = -alpha*s + beta*c;
end