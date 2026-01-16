%% =====================================================================
%% VERIFICATION - VALUE SUBSTITUTION METHOD
%% =====================================================================
% Purpose: Verify that temporal variables (Fg, Lx, Lp) reach stabilization
% Method: Check if |X(t) - X(t-1)| < tolerance for last few time steps
% Model: Fatigue simulation for military cadets (Palapes/Suksis)
% Scenario: Extreme conditions (maximum stress, zero recovery)
%% =====================================================================

clc;
clear;
close all;

%% ===================== TIME & SIM SETUP =====================
T  = 50;           % total duration (arbitrary units)
dt = 0.1;          % step size (0<dt<=1 is stable for these simple trackers)
numSteps  = T/dt;  % number of steps
time = linspace(0, T, numSteps);

%% ===================== INPUTS (EXTREME SCENARIO) =====================
% Initialize arrays to store the variables over time (all normalized [0,1])
Te = zeros(1, numSteps); % Thermal environment (heat)
Ti = zeros(1, numSteps); % Training intensity
Tc = zeros(1, numSteps); % Training complexity
Tp = zeros(1, numSteps); % Training pressure
Fi = zeros(1, numSteps); % Fluid intake
Fl = zeros(1, numSteps); % Fitness level
Hr = zeros(1, numSteps); % Heart-rate strain (0=low heart rate, 1=high heart rate)
Sh = zeros(1, numSteps); % Sleep hours (normalized)
St = zeros(1, numSteps); % Sleep timing / alignment (normalized)

% Initial values - EXTREME SCENARIO
for t = 1:numSteps
    Te(t) = 1;     % Maximum thermal stress
    Ti(t) = 1;     % Maximum training intensity
    Tc(t) = 1;     % Maximum training complexity
    Tp(t) = 1;     % Maximum training pressure
    Hr(t) = 1;     % Maximum heart rate
    Fi(t) = 0;     % No fluid intake
    Fl(t) = 0;     % No fitness
    Sh(t) = 0;     % No sleep
    St(t) = 0;     % Poor sleep timing
end

%% ===================== PARAMETERS =====================
% Task Demands
w_Td1 = 0.4;
w_Td2 = 0.6;
beta_Td = 0.7;

% Sleep Quality
delta_Sq = 0.7;

% Hydration Level
beta_Hl = 0.5;
w_Hl1 = 0.6;
w_Hl2 = 0.4;

% Physical Preparedness
beta_Pp  = 0.7;
w_Pp1 = 0.25;
w_Pp2 = 0.30;
w_Pp3 = 0.45;

% Effort Motivation
alpha_Em = 0.4;
eta_Em = 0.5;

% Generated Effort
w_Ge1 = 0.85;
w_Ge2 = 0.15;

% Recovery Effort (CALIBRATED)
w_Re1 = 0.25;
w_Re2 = 0.20;
w_Re3 = 0.20;

% Short-term Exhaustion (CALIBRATED)
w_Sx1 = 0.70;
w_Sx2 = 1.10;
gamma_Sx = 0.50;

% Short-term Experienced Pressure
w_Sp1 = 0.6;
w_Sp2 = 0.4;
beta_Sp = 0.7;

% Recovery Time
w_Rt1 = 0.5;
w_Rt2 = 0.5;
beta_Rt1 = 0.3;
beta_Rt2 = 0.3;
beta_Rt3 = 0.4;

%% ===================== STATE ARRAYS =====================
Td = zeros(1,numSteps);     % Task Demands
Hl = zeros(1,numSteps);     % Hydration Level
Sq = zeros(1,numSteps);     % Sleep Quality
Pp = zeros(1,numSteps);     % Physical Preparedness
Em = zeros(1,numSteps);     % Effort Motivation
Ge = zeros(1,numSteps);     % Generated Effort
Cp = zeros(1,numSteps);     % Critical Point
Sp = zeros(1,numSteps);     % Short-term Experienced Pressure
Sx = zeros(1,numSteps);     % Short-term Exhaustion
Re = zeros(1,numSteps);     % Recovery Effort
Lp = zeros(1,numSteps);     % Long-term Experienced Pressure
Lx = zeros(1,numSteps);     % Long-term Exhaustion
Fg = zeros(1,numSteps);     % Fatigue
Rt = zeros(1,numSteps);     % Recovery Time

%% ===================== STEP 1 (t=1) INSTANTANEOUS COMPUTE =====================
% Task Demands
Td(1) = beta_Td * Ti(1) * (w_Td1*Tc(1) + w_Td2*Tp(1)) + (1 - beta_Td) * (Te(1) * Ti(1));

% Sleep Quality
Sq(1) = delta_Sq * Sh(1) + (1-delta_Sq) * St(1);

% Hydration Level
Hl(1) = (beta_Hl * Fi(1)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(1) + w_Hl2 * Td(1)));

% Physical Preparedness
Pp(1) = beta_Pp * (w_Pp1*Hl(1) + w_Pp2*Fl(1) + w_Pp3*Sq(1)) + (1-beta_Pp) * (1 - Hr(1));

% CALIBRATED INITIAL CONDITIONS (personalized based on preparedness)
initial_stress = 1 - Pp(1);
Fg(1) = 0.35 + 0.35 * initial_stress;
Lx(1) = Fg(1);
Lp(1) = 0.40 + 0.25 * initial_stress;
Sx(1) = 0.40 + 0.30 * initial_stress;

% Critical Point
Cp(1) = Pp(1) * (1 - Sx(1));

% Effort Motivation
Em(1) = exp( -eta_Em * ((Td(1) - (1-Lp(1)) ) / alpha_Em) ^ 2 );

% Generated Effort
Ge(1) = Em(1) * (w_Ge1*Td(1) + w_Ge2*Cp(1));

% Recovery Effort
Re(1) = w_Re1 * min(1, max(0, Cp(1) - Ge(1))) + w_Re2 * Ge(1) + w_Re3 * min(1, max(0, (Pp(1) - Cp(1)) / max(Pp(1), 1e-9))) ;

% Short-term Exhaustion
Sx(1) = (w_Sx1*Fg(1) + w_Sx2*max(0, Ge(1)-Cp(1))) * (1 - gamma_Sx*Re(1));

% Short-term Experienced Pressure
Sp(1) = beta_Sp * (w_Sp1 * Sx(1) + w_Sp2 * max(0, Ge(1)-Cp(1))) + (1-beta_Sp) * (1-Pp(1));

% Recovery Time
Rt(1) = (w_Rt1*Sx(1) + w_Rt2*Fg(1)) * (1 - (beta_Rt1*Pp(1) + beta_Rt2*Sq(1) + beta_Rt3*Re(1)));

%% ===================== MAIN LOOP =====================
for t = 2:numSteps
    % -------- Instantaneous layer --------

    % Task Demands
    Td(t) = beta_Td * Ti(t) * (w_Td1*Tc(t) + w_Td2*Tp(t)) + (1 - beta_Td) * (Te(t) * Ti(t));

    % Sleep Quality
    Sq(t) = delta_Sq * Sh(t) + (1-delta_Sq) * St(t-1);

    % Hydration Level
    Hl(t) = (beta_Hl * Fi(t)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(t) + w_Hl2 * Td(t)));

    % Physical Preparedness
    Pp(t) = beta_Pp * (w_Pp1*Hl(t) + w_Pp2*Fl(t) + w_Pp3*Sq(t)) + (1-beta_Pp) * (1 - Hr(t));

    % Critical Point
    Cp(t) = Pp(t) * (1 - Sx(t-1));

    % Effort Motivation
    Em(t) = exp( -eta_Em * ((Td(t) - (1-Lp(t-1)) ) / alpha_Em) ^ 2 );

    % Generated Effort
    Ge(t) = Em(t) * (w_Ge1*Td(t) + w_Ge2*Cp(t));

    % Recovery Effort
    Re(t) = w_Re1 * min(1, max(0, Cp(t) - Ge(t))) + w_Re2 * Ge(t) + w_Re3 * min(1, max(0, (Pp(t) - Cp(t)) / max(Pp(t), 1e-9)));

    % Short-term exhaustion
    Sx(t) = (w_Sx1*Fg(t-1) + w_Sx2*max(0, Ge(t)-Cp(t))) * (1 - gamma_Sx*Re(t));

    % Short-term pressure
    Sp(t) = beta_Sp * (w_Sp1 * Sx(t-1) + w_Sp2 * max(0, Ge(t)-Cp(t))) + (1-beta_Sp) * (1-Pp(t));

    % -------- Temporal layer (first-order trackers) --------
    % Long-term experienced pressure
    Lp(t) = Lp(t-1) + (Sp(t) - Lp(t-1)) * Lp(t-1) * (1-Lp(t-1)) * dt;

    % Long-term exhaustion
    Lx(t) = Lx(t-1) + (Sx(t-1) - Lx(t-1)) * Lx(t-1) * (1-Lx(t-1)) * dt;

    % Fatigue
    Fg(t) = Fg(t-1) + (Lx(t-1) - Fg(t-1)) * Fg(t-1) * (1-Fg(t-1)) * dt;

    % Recovery Time
    Rt(t) = (w_Rt1*Sx(t-1) + w_Rt2*Fg(t-1)) * (1 - (beta_Rt1*Pp(t) + beta_Rt2*Sq(t) + beta_Rt3*Re(t)));
end

%% =====================================================================
%% VERIFICATION SECTION - VALUE SUBSTITUTION
%% =====================================================================

%% ========== CALCULATE CHANGES (DERIVATIVES) ==========
dFg = diff(Fg);  % Fg(t) - Fg(t-1)
dLx = diff(Lx);  % Lx(t) - Lx(t-1)
dLp = diff(Lp);  % Lp(t) - Lp(t-1)

%% ========== STABILIZATION CHECK ==========
tolerance = 1e-3;  % 0.001

% Determine which steps are stabilized for each variable
stable_Fg = abs(dFg) < tolerance;
stable_Lx = abs(dLx) < tolerance;
stable_Lp = abs(dLp) < tolerance;

%% ========== CONSOLE OUTPUT - STABILIZATION TABLE ==========
fprintf('\n=============================================================\n');
fprintf('         VALUE SUBSTITUTION - STABILIZATION CHECK\n');
fprintf('=============================================================\n');
fprintf('Model: Fatigue Simulation for Military Cadets\n');
fprintf('Scenario: Extreme (Te=1, Ti=1, Tc=1, Tp=1, Hr=1, Sh=0, St=0, Fi=0, Fl=0)\n');
fprintf('Time parameters: T=%d, dt=%.1f, numSteps=%d\n', T, dt, numSteps);
fprintf('Tolerance for stabilization: %.3e\n\n', tolerance);

% Select last 10 consecutive steps
last_indices = (numSteps-9):numSteps;  % e.g., [491, 492, ..., 500]
last_times = time(last_indices);  % Corresponding time values

fprintf('Last 10 consecutive time steps:\n');
fprintf('------------------------------------------------------------------------------------------------------------------------\n');
fprintf('Time      Fg          dFg         Status        Lx          dLx         Status        Lp          dLp         Status\n');
fprintf('------------------------------------------------------------------------------------------------------------------------\n');

for k = 1:length(last_indices)
    i = last_indices(k);
    idx_diff = i - 1;  % Index for diff arrays

    % Status strings
    if stable_Fg(idx_diff)
        status_Fg = 'Stabilized';
    else
        status_Fg = 'Not Stab  ';
    end

    if stable_Lx(idx_diff)
        status_Lx = 'Stabilized';
    else
        status_Lx = 'Not Stab  ';
    end

    if stable_Lp(idx_diff)
        status_Lp = 'Stabilized';
    else
        status_Lp = 'Not Stab  ';
    end

    fprintf('%-6.1f    %.6f    %.6f    %s    %.6f    %.6f    %s    %.6f    %.6f    %s\n', ...
        last_times(k), Fg(i), dFg(idx_diff), status_Fg, ...
        Lx(i), dLx(idx_diff), status_Lx, ...
        Lp(i), dLp(idx_diff), status_Lp);
end

fprintf('------------------------------------------------------------------------------------------------------------------------\n\n');

% Final summary
fprintf('CONCLUSION:\n');

% Count how many of last 4 integer time steps are stabilized
last_4_indices = last_indices(end-3:end) - 1;  % Last 4 diff indices
last_4_Fg = sum(stable_Fg(last_4_indices));
last_4_Lx = sum(stable_Lx(last_4_indices));
last_4_Lp = sum(stable_Lp(last_4_indices));

if last_4_Fg == 4
    fprintf('Fg (Fatigue):                    Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf('Fg (Fatigue):                    Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Fg, tolerance);
end

if last_4_Lx == 4
    fprintf('Lx (Long-term Exhaustion):       Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf('Lx (Long-term Exhaustion):       Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Lx, tolerance);
end

if last_4_Lp == 4
    fprintf('Lp (Long-term Pressure):         Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf('Lp (Long-term Pressure):         Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Lp, tolerance);
end

fprintf('\nFinal values:\n');
fprintf('  Fg(%d) = %.6f\n', T, Fg(end));
fprintf('  Lx(%d) = %.6f\n', T, Lx(end));
fprintf('  Lp(%d) = %.6f\n', T, Lp(end));
fprintf('=============================================================\n\n');

%% ========== FIGURE 1: FATIGUE & LONG-TERM EXHAUSTION ==========
figure('Name', 'Fatigue and Long-term Exhaustion');
plot(time, Fg, 'r-', 'LineWidth', 2); hold on;
plot(time, Lx, 'r--', 'LineWidth', 2);
xlabel('Time'); ylabel('Level');
title('Fatigue and Long-term Exhaustion');
legend('Fg (Fatigue)', 'Lx (Long-term Exhaustion)', 'Location', 'best');
grid on;
ylim([0 1.05]);

%% ========== FIGURE 2: LONG-TERM EXHAUSTION & SHORT-TERM EXHAUSTION ==========
figure('Name', 'Long-term and Short-term Exhaustion');
plot(time, Lx, 'r-', 'LineWidth', 2); hold on;
plot(time, Sx, 'r--', 'LineWidth', 2);
xlabel('Time'); ylabel('Level');
title('Long-term Exhaustion and Short-term Exhaustion');
legend('Lx (Long-term Exhaustion)', 'Sx (Short-term Exhaustion)', 'Location', 'best');
grid on;
ylim([0 1.05]);

%% ========== FIGURE 3: LONG-TERM PRESSURE & SHORT-TERM PRESSURE ==========
figure('Name', 'Long-term and Short-term Pressure');
plot(time, Lp, 'r-', 'LineWidth', 2); hold on;
plot(time, Sp, 'r--', 'LineWidth', 2);
xlabel('Time'); ylabel('Level');
title('Long-term Pressure and Short-term Pressure');
legend('Lp (Long-term Pressure)', 'Sp (Short-term Pressure)', 'Location', 'best');
grid on;
ylim([0 1.05]);

%% ========== FIGURE 4: CHANGES IN FATIGUE ==========
figure('Name', 'Changes in Fatigue');
plot(time(2:end), dFg, 'r-*', 'LineWidth', 1.5);
xlabel('Time'); ylabel('\DeltaFg');
title('Changes in Fatigue');
grid on;

%% ========== FIGURE 5: CHANGES IN LONG-TERM EXHAUSTION ==========
figure('Name', 'Changes in Long-term Exhaustion');
plot(time(2:end), dLx, 'r-*', 'LineWidth', 1.5);
xlabel('Time'); ylabel('\DeltaLx');
title('Changes in Long-term Exhaustion');
grid on;

%% ========== FIGURE 6: CHANGES IN LONG-TERM PRESSURE ==========
figure('Name', 'Changes in Long-term Pressure');
plot(time(2:end), dLp, 'r-*', 'LineWidth', 1.5);
xlabel('Time'); ylabel('\DeltaLp');
title('Changes in Long-term Pressure');
grid on;

%% ========== SAVE RESULTS TO TEXT FILE ==========
fid = fopen('verification_results.txt', 'w');

fprintf(fid, '=============================================================\n');
fprintf(fid, '         VALUE SUBSTITUTION - STABILIZATION CHECK\n');
fprintf(fid, '=============================================================\n');
fprintf(fid, 'Model: Fatigue Simulation for Military Cadets\n');
fprintf(fid, 'Scenario: Extreme (Te=1, Ti=1, Tc=1, Tp=1, Hr=1, Sh=0, St=0, Fi=0, Fl=0)\n');
fprintf(fid, 'Time parameters: T=%d, dt=%.1f, numSteps=%d\n', T, dt, numSteps);
fprintf(fid, 'Tolerance: %.3e\n\n', tolerance);

fprintf(fid, 'FULL TIME SERIES DATA:\n');
fprintf(fid, '----------------------------------------------------------------------------------------------------\n');
fprintf(fid, 'Time      Fg          dFg         Lx          dLx         Lp          dLp\n');
fprintf(fid, '----------------------------------------------------------------------------------------------------\n');
fprintf(fid, '%.1f     %.6f    -           %.6f    -           %.6f    -\n', time(1), Fg(1), Lx(1), Lp(1));
for i = 2:numSteps
    fprintf(fid, '%.1f     %.6f    %.6f    %.6f    %.6f    %.6f    %.6f\n', ...
        time(i), Fg(i), dFg(i-1), Lx(i), dLx(i-1), Lp(i), dLp(i-1));
end

fprintf(fid, '\n\nSTABILIZATION ANALYSIS (Last 10 consecutive time steps):\n');
fprintf(fid, '----------------------------------------------------------------------------------------------------\n');
fprintf(fid, 'Time      Fg          dFg         Status      Lx          dLx         Status      Lp          dLp         Status\n');
fprintf(fid, '----------------------------------------------------------------------------------------------------\n');

for k = 1:length(last_indices)
    i = last_indices(k);
    idx_diff = i - 1;

    if stable_Fg(idx_diff)
        status_Fg = 'Stabilized';
    else
        status_Fg = 'Not Stab  ';
    end

    if stable_Lx(idx_diff)
        status_Lx = 'Stabilized';
    else
        status_Lx = 'Not Stab  ';
    end

    if stable_Lp(idx_diff)
        status_Lp = 'Stabilized';
    else
        status_Lp = 'Not Stab  ';
    end

    fprintf(fid, '%-6.1f    %.6f    %.6f    %s    %.6f    %.6f    %s    %.6f    %.6f    %s\n', ...
        last_times(k), Fg(i), dFg(idx_diff), status_Fg, ...
        Lx(i), dLx(idx_diff), status_Lx, ...
        Lp(i), dLp(idx_diff), status_Lp);
end

fprintf(fid, '\n\nCONCLUSION:\n');
if last_4_Fg == 4
    fprintf(fid, 'Fg (Fatigue): Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf(fid, 'Fg (Fatigue): Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Fg, tolerance);
end

if last_4_Lx == 4
    fprintf(fid, 'Lx (Long-term Exhaustion): Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf(fid, 'Lx (Long-term Exhaustion): Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Lx, tolerance);
end

if last_4_Lp == 4
    fprintf(fid, 'Lp (Long-term Pressure): Stabilized (all last 4 steps < %.3e)\n', tolerance);
else
    fprintf(fid, 'Lp (Long-term Pressure): Not fully stabilized (%d/4 steps < %.3e)\n', last_4_Lp, tolerance);
end

fprintf(fid, '\nFinal values at t=%d:\n', T);
fprintf(fid, '  Fg = %.6f\n', Fg(end));
fprintf(fid, '  Lx = %.6f\n', Lx(end));
fprintf(fid, '  Lp = %.6f\n', Lp(end));
fprintf(fid, '=============================================================\n');

fclose(fid);

fprintf('\nResults saved to verification_results.txt\n');
fprintf('All 6 figures generated successfully.\n\n');
