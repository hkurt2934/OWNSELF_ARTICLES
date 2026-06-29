%% Sensorless IM Drive with Inverter OCF Fault Injection
clear; clc; close all;

%% ================== Simulation Settings ==================
Ts     = 1e-4;      % Sampling time [s]
Tend   = 3;         % Total simulation time [s]
tfault = 1;         % Fault injection time [s]
time   = 0:Ts:Tend;
N      = length(time);

%% ================== Motor Parameters ==================
P  = 4;             % Number of poles
Rs = 0.435;         % Stator resistance [Ohm]
Rr = 0.816;         % Rotor resistance [Ohm]
Ls = 0.002;         % Stator inductance [H]
Lr = 0.002;         % Rotor inductance [H]
Lm = 0.0693;        % Mutual inductance [H]
J  = 0.089;         % Inertia [kg.m^2]
B  = 0.001;         % Friction coefficient

Vdc = 600;          % DC link voltage [V]

%% ================== References ==================
w_ref = 1500 * 2*pi/60;   % Mechanical speed reference [rad/s]
id_ref = 0;              % Flux current reference
iq_ref = 100;            % Torque current reference

%% ================== Variables ==================
wm  = zeros(1,N);        % Mechanical speed
Te  = zeros(1,N);        % Electromagnetic torque

ia = zeros(1,N); ib = ia; ic = ia;
va = ia; vb = ia; vc = ia;

%% ================== Inverter Gate Signals ==================
Sa = ones(1,N);   % Phase-A upper switch
Sb = ones(1,N);   % Phase-B upper switch
Sc = ones(1,N);   % Phase-C upper switch

%% ================== Fault Injection ==================
for k = 1:N

    t = time(k);

    % -------- Inverter Open-Circuit Fault --------
    if t >= tfault
        Sa(k) = 0;     % Open-circuit fault in phase-A upper switch
        % For leg fault, also disable lower switch in Simulink
    end

    % -------- Phase Voltages (Simplified SVPWM) --------
    va(k) = (2*Sa(k)-1) * Vdc/2;
    vb(k) = (2*Sb(k)-1) * Vdc/2;
    vc(k) = (2*Sc(k)-1) * Vdc/2;

    % -------- Current & Torque (Simplified Placeholder) --------
    ia(k) = va(k)/Rs;
    ib(k) = vb(k)/Rs;
    ic(k) = vc(k)/Rs;

    Te(k) = 1.5 * (P/2) * Lm * iq_ref;

    % -------- Mechanical Equation --------
    if k > 1
        wm(k) = wm(k-1) + Ts*(Te(k) - B*wm(k-1))/J;
    end
end

%% ================== Plots ==================
figure;
plot(time, wm, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('Speed [rad/s]');
title('Rotor Speed Response with Inverter Open-Circuit Fault');
xline(tfault, '--r', 'Fault');

figure;
plot(time, ia, 'LineWidth', 1.5); hold on;
plot(time, ib, 'LineWidth', 1.5);
plot(time, ic, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('Phase Currents [A]');
title('Stator Phase Currents under Inverter OCF');
legend('i_a','i_b','i_c');
xline(tfault, '--r', 'Fault');
