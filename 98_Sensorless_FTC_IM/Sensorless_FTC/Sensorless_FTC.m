%% Induction Motor + FOC + FTC Simulation
clear; clc; close all;

%% Motor Parameters (30 kW EV Induction Motor)
P = 4;                % Pole pairs (2 pole-pairs = 4 poles)
Rs = 0.435;           % Stator resistance [Ohm]
Rr = 0.816;           % Rotor resistance [Ohm]
Ls = 0.002;           % Stator inductance [H]
Lr = 0.002;           % Rotor inductance [H]
Lm = 0.0693;          % Mutual inductance [H]
J  = 0.089;           % Inertia [kg.m^2]
B  = 0.001;           % Friction coeff

Vdc = 600;            % DC link voltage [V]
fsw = 10000;          % Switching frequency [Hz]

%% Simulation Setup
Ts  = 1e-5;           % Simulation step
Tend = 1;             % End time [s]
time = 0:Ts:Tend;

% References
w_ref = 1500 * 2*pi/60;   % 1500 rpm -> rad/s
iq_ref = 100;             % Torque-producing current [A]
id_ref = 0;               % Flux-producing current [A]

%% Controller Gains (tuned roughly)
Kp_speed = 1; Ki_speed = 50;
Kp_curr = 5;  Ki_curr = 200;

%% Variables
n = length(time);
ids = zeros(1,n); iqs = zeros(1,n);
idr = zeros(1,n); iqr = zeros(1,n);
theta_r = zeros(1,n); omega_r = zeros(1,n);
Te = zeros(1,n);

iqs_ref = zeros(1,n); id_ref_vec = id_ref*ones(1,n);

fault_flag = false;

%% MRAS Observer variables
omega_hat = zeros(1,n);   % estimated speed
K_adapt = 500;

%% Loop
for k = 2:n
    
    % --- Fault injection at t=0.5s ---
    if time(k) > 0.5 && ~fault_flag
        disp('>>> Fault inserted: S1 open at t=0.5s');
        fault_flag = true;
    end
    
    % --- Control ---
    % Speed controller (outer loop) -> iqs_ref
    e_speed = w_ref - omega_r(k-1);
    iqs_ref(k) = iqs_ref(k-1) + Ts*(Kp_speed*e_speed + Ki_speed*e_speed);
    
    % Current controllers (inner loop) -> simple proportional (dq control)
    vds = Kp_curr*(id_ref - ids(k-1));
    vqs = Kp_curr*(iqs_ref(k) - iqs(k-1));
    
    % --- Fault effect ---
    if fault_flag
        % simulate open S1 (phase A upper) -> reduce vds
        vds = 0.7*vds; % degraded voltage
    end
    
    % --- Motor model (Euler integration) ---
    % Flux linkages
    lambda_ds = Ls*ids(k-1) + Lm*idr(k-1);
    lambda_qs = Ls*iqs(k-1) + Lm*iqr(k-1);
    lambda_dr = Lr*idr(k-1) + Lm*ids(k-1);
    lambda_qr = Lr*iqr(k-1) + Lm*iqs(k-1);
    
    % Stator voltage equations (simplified)
    dids = (vds - Rs*ids(k-1) + omega_r(k-1)*lambda_qs) / Ls;
    diqs = (vqs - Rs*iqs(k-1) - omega_r(k-1)*lambda_ds) / Ls;
    
    % Rotor dynamics
    didr = (-Rr*idr(k-1) + (omega_r(k-1)-omega_hat(k-1))*lambda_qr) / Lr;
    diqr = (-Rr*iqr(k-1) - (omega_r(k-1)-omega_hat(k-1))*lambda_dr) / Lr;
    
    % Torque
    Te(k) = (3/2)*(P/2)*Lm*(ids(k-1)*iqr(k-1) - iqs(k-1)*idr(k-1));
    
    % Mechanical equation
    domega = (Te(k) - B*omega_r(k-1))/J;
    
    % Integrate
    ids(k) = ids(k-1) + Ts*dids;
    iqs(k) = iqs(k-1) + Ts*diqs;
    idr(k) = idr(k-1) + Ts*didr;
    iqr(k) = iqr(k-1) + Ts*diqr;
    omega_r(k) = omega_r(k-1) + Ts*domega;
    theta_r(k) = theta_r(k-1) + Ts*omega_r(k);
    
    % --- MRAS speed estimation ---
    eps_mras = (lambda_ds - (Ls*ids(k)+Lm*idr(k))) + ...
               (lambda_qs - (Ls*iqs(k)+Lm*iqr(k)));
    omega_hat(k) = omega_hat(k-1) + Ts*K_adapt*eps_mras;
end

%% Post Processing Metrics
% Torque ripple
torque_ripple = (max(Te) - min(Te))/mean(Te)*100;

% Speed error
speed_err = mean(abs(w_ref - omega_r));

fprintf('Torque ripple: %.2f %%\n', torque_ripple);
fprintf('Mean speed error: %.2f rad/s\n', speed_err);

%% Plots
figure;
subplot(3,1,1); plot(time,omega_r*60/(2*pi),'b','LineWidth',1.2); hold on;
yline(w_ref*60/(2*pi),'r--'); ylabel('Speed [rpm]'); grid on;

subplot(3,1,2); plot(time,Te,'k','LineWidth',1.2);
ylabel('Torque [Nm]'); grid on;

subplot(3,1,3); plot(time,ids,time,iqs); legend('i_d','i_q');
xlabel('Time [s]'); ylabel('Currents [A]'); grid on;
