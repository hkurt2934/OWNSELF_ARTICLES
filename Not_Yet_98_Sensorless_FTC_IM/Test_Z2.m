%% ========================================================================
%  SENSORLESS FOC INDUCTION MOTOR WITH INVERTER LEG FAULT (NO FTC)
%  ========================================================================
%  Description:
%    Simulation of a three-phase induction motor drive using indirect
%    rotor flux-oriented control (IRFOC) with MRAS speed observer under
%    a single inverter leg open-circuit fault in Phase A.
%
%    FAULT CHARACTERISTICS:
%    - Both switches in Phase A leg are open (cannot conduct)
%    - Phase A current forced to zero after fault
%    - No fault-tolerant control or reconfiguration
%    - Controllers continue operating without modification
%
%  Author: Motor Control Research Group
%  Date: 2025
%  Purpose: Fault case for comparison with healthy operation and FTC
%  ========================================================================

clear all; close all; clc;

%% ========================================================================
%  1. SIMULATION PARAMETERS
%  ========================================================================
Ts = 100e-6;                    % Sampling time [s] (100 microseconds)
T_sim = 3;                      % Simulation time [s]
time = 0:Ts:T_sim;             % Time vector
N = length(time);              % Number of samples

% Fault Parameters
T_fault = 1.0;                 % Fault injection time [s]
fault_index = round(T_fault/Ts);  % Fault sample index

%% ========================================================================
%  2. INDUCTION MOTOR PARAMETERS (Standard 2.2 kW Motor)
%  ========================================================================
% Electrical Parameters
P = 4;                         % Number of poles
Rs = 1.115;                    % Stator resistance [Ohm]
Rr = 1.083;                    % Rotor resistance [Ohm]
Lls = 0.005974;               % Stator leakage inductance [H]
Llr = 0.005974;               % Rotor leakage inductance [H]
Lm = 0.2037;                  % Magnetizing inductance [H]
Ls = Lls + Lm;                % Stator self-inductance [H]
Lr = Llr + Lm;                % Rotor self-inductance [H]

% Time Constants
Tr = Lr/Rr;                   % Rotor time constant [s]
sigma = 1 - (Lm^2)/(Ls*Lr);   % Leakage coefficient

% Mechanical Parameters
J = 0.0131;                   % Moment of inertia [kg.m^2]
B = 0.002;                    % Friction coefficient [N.m.s/rad]
T_rated = 14.5;               % Rated torque [N.m]

%% ========================================================================
%  3. CONTROL PARAMETERS
%  ========================================================================
% Reference Values
omega_ref_rpm = 1440;         % Reference speed [rpm]
omega_ref = omega_ref_rpm * (2*pi/60);  % Reference speed [rad/s]
i_ds_ref = 5.0;              % d-axis current reference (flux producing) [A]

% Load Torque Profile
T_load = 10.0 * ones(1, N);  % Constant load torque [N.m]

% PI Controller Gains - Speed Loop
Kp_speed = 0.5;              % Proportional gain
Ki_speed = 5.0;              % Integral gain

% PI Controller Gains - Current Loops (d and q axes)
Kp_id = 10.0;                % d-axis current proportional gain
Ki_id = 100.0;               % d-axis current integral gain
Kp_iq = 10.0;                % q-axis current proportional gain
Ki_iq = 100.0;               % q-axis current integral gain

% MRAS Observer Gains
Kp_mras = 50;                % MRAS proportional gain
Ki_mras = 300;               % MRAS integral gain

%% ========================================================================
%  4. INITIALIZATION OF STATE VARIABLES
%  ========================================================================
% Motor States (in synchronous dq reference frame)
i_ds = 0;                    % d-axis stator current [A]
i_qs = 0;                    % q-axis stator current [A]
psi_dr = 0;                  % d-axis rotor flux [Wb]
psi_qr = 0;                  % q-axis rotor flux [Wb]
omega_r = 0;                 % Actual rotor electrical speed [rad/s]
omega_m = 0;                 % Actual rotor mechanical speed [rad/s]
theta_e = 0;                 % Rotor flux angle [rad]

% Three-phase currents (in abc frame)
i_as = 0;                    % Phase A current [A]
i_bs = 0;                    % Phase B current [A]
i_cs = 0;                    % Phase C current [A]

% MRAS Observer States (in stationary αβ reference frame)
psi_alpha_ref = 0;           % α-axis rotor flux - reference model [Wb]
psi_beta_ref = 0;            % β-axis rotor flux - reference model [Wb]
psi_alpha_adj = 0;           % α-axis rotor flux - adjustable model [Wb]
psi_beta_adj = 0;            % β-axis rotor flux - adjustable model [Wb]
omega_r_est = 0;             % Estimated rotor electrical speed [rad/s]
omega_m_est = 0;             % Estimated rotor mechanical speed [rad/s]

% PI Controller Integral States
int_speed = 0;               % Speed controller integral
int_id = 0;                  % d-axis current controller integral
int_iq = 0;                  % q-axis current controller integral
int_mras = 0;                % MRAS observer integral

% Fault Status Flag
fault_active = false;        % Flag indicating fault status

%% ========================================================================
%  5. STORAGE ARRAYS FOR RESULTS
%  ========================================================================
omega_m_actual = zeros(1, N);
omega_m_estimated = zeros(1, N);
Te_array = zeros(1, N);
i_ds_array = zeros(1, N);
i_qs_array = zeros(1, N);
i_as_array = zeros(1, N);
i_bs_array = zeros(1, N);
i_cs_array = zeros(1, N);
fault_status = zeros(1, N);   % Fault status indicator

%% ========================================================================
%  6. MAIN SIMULATION LOOP
%  ========================================================================
fprintf('Starting simulation with inverter leg fault...\n');
fprintf('Fault will be injected at t = %.2f seconds\n', T_fault);

for k = 1:N
    %% --------------------------------------------------------------------
    %  6.1 FAULT INJECTION LOGIC
    %  --------------------------------------------------------------------
    if k >= fault_index && ~fault_active
        fault_active = true;
        fprintf('FAULT INJECTED at t = %.3f s: Phase A inverter leg open-circuit\n', time(k));
    end
    
    fault_status(k) = fault_active;
    
    %% --------------------------------------------------------------------
    %  6.2 SPEED CONTROLLER (Outer Loop) - Unchanged
    %  --------------------------------------------------------------------
    % Speed error
    error_speed = omega_ref - omega_m_est;
    
    % PI controller for speed
    int_speed = int_speed + error_speed * Ts;
    i_qs_ref = Kp_speed * error_speed + Ki_speed * int_speed;
    
    % Limit q-axis current reference
    i_qs_ref = max(min(i_qs_ref, 15), -15);
    
    %% --------------------------------------------------------------------
    %  6.3 CURRENT CONTROLLERS (Inner Loop) - Unchanged
    %  --------------------------------------------------------------------
    % d-axis current controller
    error_id = i_ds_ref - i_ds;
    int_id = int_id + error_id * Ts;
    v_ds_ref = Kp_id * error_id + Ki_id * int_id;
    
    % q-axis current controller
    error_iq = i_qs_ref - i_qs;
    int_iq = int_iq + error_iq * Ts;
    v_qs_ref = Kp_iq * error_iq + Ki_iq * int_iq;
    
    % Voltage limiting (simple saturation)
    V_max = 400;  % Maximum voltage [V]
    v_mag = sqrt(v_ds_ref^2 + v_qs_ref^2);
    if v_mag > V_max
        v_ds_ref = v_ds_ref * (V_max / v_mag);
        v_qs_ref = v_qs_ref * (V_max / v_mag);
    end
    
    %% --------------------------------------------------------------------
    %  6.4 INDIRECT FIELD ORIENTATION - Unchanged
    %  --------------------------------------------------------------------
    % Slip frequency calculation
    if abs(psi_dr) > 0.01
        omega_sl = (Lm / (Tr * psi_dr)) * i_qs;
    else
        omega_sl = 0;
    end
    
    % Synchronous frequency
    omega_e = omega_r_est + omega_sl;
    
    % Rotor flux angle integration
    theta_e = theta_e + omega_e * Ts;
    theta_e = mod(theta_e, 2*pi);  % Keep angle in [0, 2π]
    
    %% --------------------------------------------------------------------
    %  6.5 INVERSE PARK TRANSFORMATION (dq to αβ) - Unchanged
    %  --------------------------------------------------------------------
    v_alpha = v_ds_ref * cos(theta_e) - v_qs_ref * sin(theta_e);
    v_beta = v_ds_ref * sin(theta_e) + v_qs_ref * cos(theta_e);
    
    %% --------------------------------------------------------------------
    %  6.6 INVERSE CLARKE TRANSFORMATION (αβ to abc) - Unchanged
    %  --------------------------------------------------------------------
    v_as = v_alpha;
    v_bs = -0.5 * v_alpha + (sqrt(3)/2) * v_beta;
    v_cs = -0.5 * v_alpha - (sqrt(3)/2) * v_beta;
    
    %% --------------------------------------------------------------------
    %  6.7 FAULT IMPLEMENTATION: Phase A Open-Circuit
    %  --------------------------------------------------------------------
    if fault_active
        % Phase A inverter leg is open-circuit: voltage cannot be applied
        v_as = 0;
        
        % Phase A current forced to zero (open circuit)
        i_as = 0;
        
        % Kirchhoff's current law: i_a + i_b + i_c = 0
        % Since i_a = 0, we have: i_b + i_c = 0
        % The motor operates with two phases only
    end
    
    %% --------------------------------------------------------------------
    %  6.8 CLARKE TRANSFORMATION (abc to αβ) - Unchanged Equations
    %  --------------------------------------------------------------------
    % Note: During fault, i_as = 0, so the transformation receives [0, i_bs, i_cs]
    % The Clarke transformation is NOT reconfigured for fault tolerance
    i_alpha = i_as;
    i_beta = (1/sqrt(3)) * (i_bs - i_cs);
    
    %% --------------------------------------------------------------------
    %  6.9 PARK TRANSFORMATION (αβ to dq) - Unchanged
    %  --------------------------------------------------------------------
    i_ds = i_alpha * cos(theta_e) + i_beta * sin(theta_e);
    i_qs = -i_alpha * sin(theta_e) + i_beta * cos(theta_e);
    
    %% --------------------------------------------------------------------
    %  6.10 MOTOR MODEL - STATE DERIVATIVES
    %  --------------------------------------------------------------------
    % Transform applied voltages to αβ frame
    v_alpha_applied = v_as;
    v_beta_applied = (1/sqrt(3)) * (v_bs - v_cs);
    
    % Transform to dq frame
    v_ds_applied = v_alpha_applied * cos(theta_e) + v_beta_applied * sin(theta_e);
    v_qs_applied = -v_alpha_applied * sin(theta_e) + v_beta_applied * cos(theta_e);
    
    % Current derivatives (in dq frame)
    di_ds = (1/(sigma*Ls)) * (v_ds_applied - Rs*i_ds + omega_e*sigma*Ls*i_qs ...
            + (Lm/Lr)*(Rr*psi_dr + omega_e*psi_qr - (Lm/Tr)*i_ds));
    
    di_qs = (1/(sigma*Ls)) * (v_qs_applied - Rs*i_qs - omega_e*sigma*Ls*i_ds ...
            + (Lm/Lr)*(Rr*psi_qr - omega_e*psi_dr - (Lm/Tr)*i_qs));
    
    % Rotor flux derivatives (in dq frame)
    dpsi_dr = (Lm/Tr)*i_ds - (1/Tr)*psi_dr + (omega_e - omega_r)*psi_qr;
    dpsi_qr = (Lm/Tr)*i_qs - (1/Tr)*psi_qr - (omega_e - omega_r)*psi_dr;
    
    % Electromagnetic torque
    Te = (3/2) * (P/2) * Lm * (i_qs * psi_dr - i_ds * psi_qr);
    
    % Mechanical speed derivative
    domega_m = (1/J) * (Te - T_load(k) - B*omega_m);
    omega_r = (P/2) * omega_m;  % Electrical speed
    
    %% --------------------------------------------------------------------
    %  6.11 MRAS SPEED OBSERVER - Unchanged Algorithm
    %  --------------------------------------------------------------------
    % Reference Model (Voltage Model)
    dpsi_alpha_ref = (Lm/Tr) * (v_alpha_applied - Rs*i_alpha) - (1/Tr) * psi_alpha_ref;
    dpsi_beta_ref = (Lm/Tr) * (v_beta_applied - Rs*i_beta) - (1/Tr) * psi_beta_ref;
    
    % Adjustable Model (Current Model with estimated speed)
    dpsi_alpha_adj = -(1/Tr) * psi_alpha_adj + (Lm/Tr) * i_alpha ...
                     - omega_r_est * psi_beta_adj;
    dpsi_beta_adj = -(1/Tr) * psi_beta_adj + (Lm/Tr) * i_beta ...
                    + omega_r_est * psi_alpha_adj;
    
    % Adaptation Mechanism (Cross-product error)
    epsilon = psi_alpha_ref * psi_beta_adj - psi_beta_ref * psi_alpha_adj;
    
    % PI controller for speed estimation
    int_mras = int_mras + epsilon * Ts;
    omega_r_est_new = Kp_mras * epsilon + Ki_mras * int_mras;
    
    % Low-pass filter on estimated speed (for stability)
    tau_filter = 0.01;  % Filter time constant
    omega_r_est = omega_r_est + (omega_r_est_new - omega_r_est) * (Ts/tau_filter);
    omega_m_est = omega_r_est / (P/2);
    
    %% --------------------------------------------------------------------
    %  6.12 STATE INTEGRATION (Euler Method)
    %  --------------------------------------------------------------------
    i_ds = i_ds + di_ds * Ts;
    i_qs = i_qs + di_qs * Ts;
    psi_dr = psi_dr + dpsi_dr * Ts;
    psi_qr = psi_qr + dpsi_qr * Ts;
    omega_m = omega_m + domega_m * Ts;
    
    psi_alpha_ref = psi_alpha_ref + dpsi_alpha_ref * Ts;
    psi_beta_ref = psi_beta_ref + dpsi_beta_ref * Ts;
    psi_alpha_adj = psi_alpha_adj + dpsi_alpha_adj * Ts;
    psi_beta_adj = psi_beta_adj + dpsi_beta_adj * Ts;
    
    %% --------------------------------------------------------------------
    %  6.13 RECONSTRUCT THREE-PHASE CURRENTS FROM dq
    %  --------------------------------------------------------------------
    % Inverse Park (dq to αβ)
    i_alpha_reconstructed = i_ds * cos(theta_e) - i_qs * sin(theta_e);
    i_beta_reconstructed = i_ds * sin(theta_e) + i_qs * cos(theta_e);
    
    % Inverse Clarke (αβ to abc)
    if ~fault_active
        % Normal operation: standard inverse Clarke
        i_as = i_alpha_reconstructed;
        i_bs = -0.5 * i_alpha_reconstructed + (sqrt(3)/2) * i_beta_reconstructed;
        i_cs = -0.5 * i_alpha_reconstructed - (sqrt(3)/2) * i_beta_reconstructed;
    else
        % During fault: Phase A is forced to zero
        i_as = 0;
        % With i_a = 0 and constraint i_b + i_c = 0:
        i_bs = -0.5 * i_alpha_reconstructed + (sqrt(3)/2) * i_beta_reconstructed;
        i_cs = -i_bs;  % Enforce i_b + i_c = 0
    end
    
    %% --------------------------------------------------------------------
    %  6.14 STORE RESULTS
    %  --------------------------------------------------------------------
    omega_m_actual(k) = omega_m * (60/(2*pi));      % Convert to rpm
    omega_m_estimated(k) = omega_m_est * (60/(2*pi)); % Convert to rpm
    Te_array(k) = Te;
    i_ds_array(k) = i_ds;
    i_qs_array(k) = i_qs;
    i_as_array(k) = i_as;
    i_bs_array(k) = i_bs;
    i_cs_array(k) = i_cs;
end

fprintf('Simulation completed successfully.\n');

%% ========================================================================
%  7. RESULTS VISUALIZATION
%  ========================================================================
fprintf('Generating plots...\n');

figure('Position', [100, 100, 1200, 900], 'Color', 'w');

% Plot 1: Rotor Speed (Actual vs Estimated)
subplot(4,1,1);
plot(time, omega_m_actual, 'b-', 'LineWidth', 1.5); hold on;
plot(time, omega_m_estimated, 'r--', 'LineWidth', 1.5);
plot(time, omega_ref_rpm*ones(1,N), 'k:', 'LineWidth', 1);
% Mark fault injection time
xline(T_fault, 'r--', 'LineWidth', 2, 'Label', 'Fault Injection', ...
      'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
grid on;
xlabel('Time [s]', 'FontSize', 11);
ylabel('Speed [rpm]', 'FontSize', 11);
title('Rotor Mechanical Speed: Actual vs Estimated (With Inverter Fault)', ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Actual Speed', 'Estimated Speed (MRAS)', 'Reference Speed', ...
       'Location', 'northeast');

% Plot 2: Electromagnetic Torque
subplot(4,1,2);
plot(time, Te_array, 'b-', 'LineWidth', 1.5); hold on;
plot(time, T_load, 'r--', 'LineWidth', 1.5);
xline(T_fault, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 11);
ylabel('Torque [N.m]', 'FontSize', 11);
title('Electromagnetic Torque and Load Torque (With Inverter Fault)', ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Electromagnetic Torque', 'Load Torque', 'Location', 'northeast');

% Plot 3: dq-axis Stator Currents
subplot(4,1,3);
plot(time, i_ds_array, 'b-', 'LineWidth', 1.5); hold on;
plot(time, i_qs_array, 'r-', 'LineWidth', 1.5);
xline(T_fault, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 11);
ylabel('Current [A]', 'FontSize', 11);
title('Stator Currents in dq Reference Frame', ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('i_{ds} (Flux-producing)', 'i_{qs} (Torque-producing)', ...
       'Location', 'northeast');

% Plot 4: Three-Phase Stator Currents (abc frame)
subplot(4,1,4);
plot(time, i_as_array, 'r-', 'LineWidth', 1.5); hold on;
plot(time, i_bs_array, 'g-', 'LineWidth', 1.5);
plot(time, i_cs_array, 'b-', 'LineWidth', 1.5);
xline(T_fault, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 11);
ylabel('Current [A]', 'FontSize', 11);
title('Three-Phase Stator Currents (Phase A Fault at t=1s)', ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('i_a (Faulty)', 'i_b', 'i_c', 'Location', 'northeast');

% Main title
sgtitle('Sensorless FOC with Single Inverter Leg Fault (Phase A) - No FTC', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%  8. DETAILED ANALYSIS PLOTS
%  ========================================================================
figure('Position', [150, 150, 1200, 900], 'Color', 'w');

% Plot 5: Three-Phase Currents - Pre-Fault Zoomed
subplot(3,2,1);
t_start = 0.8; t_end = 0.95;
idx = find(time >= t_start & time <= t_end);
plot(time(idx), i_as_array(idx), 'r-', 'LineWidth', 1.2); hold on;
plot(time(idx), i_bs_array(idx), 'g-', 'LineWidth', 1.2);
plot(time(idx), i_cs_array(idx), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Current [A]', 'FontSize', 10);
title('Pre-Fault: Three-Phase Currents (Balanced)', 'FontSize', 11, 'FontWeight', 'bold');
legend('i_a', 'i_b', 'i_c', 'Location', 'northeast');

% Plot 6: Three-Phase Currents - Post-Fault Zoomed
subplot(3,2,2);
t_start = 1.5; t_end = 1.65;
idx = find(time >= t_start & time <= t_end);
plot(time(idx), i_as_array(idx), 'r-', 'LineWidth', 1.2); hold on;
plot(time(idx), i_bs_array(idx), 'g-', 'LineWidth', 1.2);
plot(time(idx), i_cs_array(idx), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Current [A]', 'FontSize', 10);
title('Post-Fault: Three-Phase Currents (Unbalanced)', 'FontSize', 11, 'FontWeight', 'bold');
legend('i_a (=0)', 'i_b', 'i_c', 'Location', 'northeast');

% Plot 7: Torque Ripple - Pre-Fault
subplot(3,2,3);
t_start = 0.8; t_end = 0.95;
idx = find(time >= t_start & time <= t_end);
plot(time(idx), Te_array(idx), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Torque [N.m]', 'FontSize', 10);
title('Pre-Fault: Electromagnetic Torque (Low Ripple)', 'FontSize', 11, 'FontWeight', 'bold');

% Plot 8: Torque Ripple - Post-Fault
subplot(3,2,4);
t_start = 1.5; t_end = 1.65;
idx = find(time >= t_start & time <= t_end);
plot(time(idx), Te_array(idx), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Torque [N.m]', 'FontSize', 10);
title('Post-Fault: Electromagnetic Torque (High Ripple)', 'FontSize', 11, 'FontWeight', 'bold');

% Plot 9: Speed Estimation Error
subplot(3,2,5);
speed_error = omega_m_actual - omega_m_estimated;
plot(time, speed_error, 'b-', 'LineWidth', 1.5); hold on;
xline(T_fault, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Speed Error [rpm]', 'FontSize', 10);
title('Speed Estimation Error (Actual - Estimated)', 'FontSize', 11, 'FontWeight', 'bold');

% Plot 10: Speed Degradation
subplot(3,2,6);
speed_degradation = omega_ref_rpm - omega_m_actual;
plot(time, speed_degradation, 'b-', 'LineWidth', 1.5); hold on;
xline(T_fault, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 10);
ylabel('Speed Drop [rpm]', 'FontSize', 10);
title('Speed Degradation from Reference', 'FontSize', 11, 'FontWeight', 'bold');

sgtitle('Detailed Fault Impact Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%  9. PERFORMANCE METRICS AND FAULT IMPACT ANALYSIS
%  ========================================================================
fprintf('\n========================================\n');
fprintf('FAULT IMPACT ANALYSIS\n');
fprintf('========================================\n');

% Pre-fault steady-state (0.5s to 1.0s)
pre_fault_idx = find(time >= 0.5 & time < T_fault);

% Post-fault steady-state (2.0s to 3.0s)
post_fault_idx = find(time >= 2.0);

% Speed Performance
speed_pre = mean(omega_m_actual(pre_fault_idx));
speed_post = mean(omega_m_actual(post_fault_idx));
speed_drop = speed_pre - speed_post;
speed_drop_percent = (speed_drop / speed_pre) * 100;

fprintf('\nSpeed Performance:\n');
fprintf('  Pre-Fault Average Speed: %.2f rpm\n', speed_pre);
fprintf('  Post-Fault Average Speed: %.2f rpm\n', speed_post);
fprintf('  Speed Drop: %.2f rpm (%.2f%%)\n', speed_drop, speed_drop_percent);

% Torque Ripple
Te_pre = Te_array(pre_fault_idx);
Te_post = Te_array(post_fault_idx);
ripple_pre = std(Te_pre);
ripple_post = std(Te_post);
ripple_increase = ((ripple_post - ripple_pre) / ripple_pre) * 100;

fprintf('\nTorque Performance:\n');
fprintf('  Pre-Fault Average Torque: %.2f N.m\n', mean(Te_pre));
fprintf('  Post-Fault Average Torque: %.2f N.m\n', mean(Te_post));
fprintf('  Pre-Fault Torque Ripple (std): %.3f N.m\n', ripple_pre);
fprintf('  Post-Fault Torque Ripple (std): %.3f N.m\n', ripple_post);
fprintf('  Ripple Increase: %.2f%%\n', ripple_increase);

% Current Unbalance
i_a_rms_post = rms(i_as_array(post_fault_idx));
i_b_rms_post = rms(i_bs_array(post_fault_idx));
i_c_rms_post = rms(i_cs_array(post_fault_idx));

fprintf('\nCurrent Analysis (Post-Fault RMS):\n');
fprintf('  Phase A Current (RMS): %.3f A (forced to zero)\n', i_a_rms_post);
fprintf('  Phase B Current (RMS): %.3f A\n', i_b_rms_post);
fprintf('  Phase C Current (RMS): %.3f A\n', i_c_rms_post);

% Speed Estimation Accuracy
est_error_pre = mean(abs(omega_m_actual(pre_fault_idx) - omega_m_estimated(pre_fault_idx)));
est_error_post = mean(abs(omega_m_actual(post_fault_idx) - omega_m_estimated(post_fault_idx)));

fprintf('\nSpeed Estimation Performance:\n');
fprintf('  Pre-Fault Estimation Error: %.2f rpm\n', est_error_pre);
fprintf('  Post-Fault Estimation Error: %.2f rpm\n', est_error_post);

fprintf('\n========================================\n');
fprintf('FAULT CHARACTERISTICS:\n');
fprintf('  Fault Type: Single inverter leg open-circuit\n');
fprintf('  Faulty Phase: Phase A\n');
fprintf('  Fault Injection Time: %.2f s\n', T_fault);
fprintf('  Fault-Tolerant Control: DISABLED\n');
fprintf('  Control Reconfiguration: NONE\n');
fprintf('========================================\n');

fprintf('\nSimulation completed. Plots generated.\n');
fprintf('This fault case can be compared with:\n');
fprintf('  1. Healthy operation (baseline)\n');
fprintf('  2. Fault-tolerant control strategies\n');

%% END OF SCRIPT