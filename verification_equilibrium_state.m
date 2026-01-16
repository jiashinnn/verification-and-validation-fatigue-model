%% =====================================================================
%% VERIFICATION: EQUILIBRIUM STATE METHOD
%% =====================================================================
% Purpose: Verify that temporal variables (Lp, Lx, Fg) can be forced to
%          equilibrium states and observe how the rest of the model responds.
%
% Method: Force Lp, Lx, Fg to specific values (0 or 1) and calculate all
%         other variables normally to check for logical consistency.
%
% Test Cases:
%   Case 1: Lp=0 ∧ Lx=0 ∧ Fg=0 (No stress/fatigue)
%   Case 2: Lp=1 ∧ Lx=1 ∧ Fg=1 (Complete saturation)
%   Case 3: Lp=Sp, Lx=Sx, Fg=Lx (Moderate Equilibrium)
%   Case 4: Lp=Sp, Lx=0, Fg=0 (Pressure Equilibrium, No Exhaustion)
%   Case 5: Lp=1, Lx=1, Fg=Lx (Maximum Saturation, Fg tracks Lx)
%
%% =====================================================================

clc;
clear;
close all;

%% ===================== TIME & SIM SETUP =====================
T  = 50;           % total duration
dt = 0.1;          % step size
numSteps  = T/dt;  % number of steps
time = linspace(0, T, numSteps);

%% ===================== INPUTS (WILL BE SET PER CASE) =====================
% Initialize arrays (all normalized [0,1])
% These will be populated for each test case with appropriate values

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

%% ===================== DEFINE TEST CASES =====================
test_cases = struct();

% Case 1: All Zero (No stress/fatigue)
test_cases(1).name = 'Case 1: Lp=0, Lx=0, Fg=0 (No Stress/Fatigue)';
test_cases(1).Lp_fixed = 0;
test_cases(1).Lx_fixed = 0;
test_cases(1).Fg_fixed = 0;

% Case 2: All Maximum (Complete saturation)
test_cases(2).name = 'Case 2: Lp=1, Lx=1, Fg=1 (Complete Saturation)';
test_cases(2).Lp_fixed = 1;
test_cases(2).Lx_fixed = 1;
test_cases(2).Fg_fixed = 1;

% Case 3: Equilibrium consistency test (Lp=Sp, Lx=Sx, Fg=Lx)
% Note: We cannot force Sp and Sx, only verify they match forced Lp and Lx
% This case tests if Fg=Lx equilibrium is maintained with moderate stress
test_cases(3).name = 'Case 3: Lp=Sp, Lx=Sx, Fg=Lx (Moderate Equilibrium)';
test_cases(3).Lp_fixed = 0.5;  % Moderate pressure
test_cases(3).Lx_fixed = 0.5;  % Moderate exhaustion
test_cases(3).Fg_fixed = 0.5;  % Fg = Lx (fatigue tracks exhaustion)

% Case 4: Lp=Sp with no exhaustion/fatigue
test_cases(4).name = 'Case 4: Lp=Sp, Lx=0, Fg=0 (Pressure Equilibrium, No Exhaustion)';
test_cases(4).Lp_fixed = 0;  
test_cases(4).Lx_fixed = 0;
test_cases(4).Fg_fixed = 0;

% Case 5: Maximum exhaustion with Fg=Lx
test_cases(5).name = 'Case 5: Lp=1, Lx=1, Fg=Lx (Maximum Saturation, Fg tracks Lx)';
test_cases(5).Lp_fixed = 1;
test_cases(5).Lx_fixed = 1;
test_cases(5).Fg_fixed = 1;  % Fg=Lx (both at maximum)




%% ===================== RUN ALL TEST CASES =====================
num_cases = length(test_cases);
results = cell(num_cases, 1);

for case_idx = 1:num_cases
    fprintf('\n========================================\n');
    fprintf('Running: %s\n', test_cases(case_idx).name);
    fprintf('========================================\n');

    % Get fixed values for this case
    Lp_fixed = test_cases(case_idx).Lp_fixed;
    Lx_fixed = test_cases(case_idx).Lx_fixed;
    Fg_fixed = test_cases(case_idx).Fg_fixed;

    %% -------- SET INPUTS ACCORDING TO CASE --------
    % Initialize input arrays
    Te = zeros(1, numSteps); % Thermal environment
    Ti = zeros(1, numSteps); % Training intensity
    Tc = zeros(1, numSteps); % Training complexity
    Tp = zeros(1, numSteps); % Training pressure
    Fi = zeros(1, numSteps); % Fluid intake
    Fl = zeros(1, numSteps); % Fitness level
    Hr = zeros(1, numSteps); % Heart-rate strain
    Sh = zeros(1, numSteps); % Sleep hours
    St = zeros(1, numSteps); % Sleep timing

    if case_idx == 1
        % Case 1: Lp=0, Lx=0, Fg=0 (No stress/fatigue)
        % Good conditions: minimal stress, maximum recovery
        fprintf('Input Scenario: Good conditions (minimal stress, maximum recovery)\n');
        for t = 1:numSteps
            Te(t) = 0;     % No thermal stress
            Ti(t) = 0;     % No training intensity
            Tc(t) = 0;     % No complexity
            Tp(t) = 0;     % No pressure
            Hr(t) = 0;     % No heart rate strain
            Fi(t) = 1;     % Maximum fluid intake
            Fl(t) = 1;     % Maximum fitness level
            Sh(t) = 1;     % Maximum sleep hours
            St(t) = 1;     % Optimal sleep timing
        end

    elseif case_idx == 2
        % Case 2: Lp=1, Lx=1, Fg=1 (Complete saturation)
        % Extreme conditions: maximum stress, zero recovery
        fprintf('Input Scenario: Extreme conditions (maximum stress, zero recovery)\n');
        for t = 1:numSteps
            Te(t) = 1;     % Maximum thermal stress
            Ti(t) = 1;     % Maximum training intensity
            Tc(t) = 1;     % Maximum complexity
            Tp(t) = 1;     % Maximum pressure
            Hr(t) = 1;     % Maximum heart rate strain
            Fi(t) = 0;     % No fluid intake
            Fl(t) = 0;     % No fitness level
            Sh(t) = 0;     % No sleep hours
            St(t) = 0;     % Poor sleep timing
        end

    elseif case_idx == 3
        % Case 3: Lp=Sp, Lx=Sx, Fg=Lx (Moderate equilibrium)
        % Balanced moderate conditions to test equilibrium consistency
        fprintf('Input Scenario: Balanced moderate conditions (testing Lp=Sp, Lx=Sx, Fg=Lx)\n');
        for t = 1:numSteps
            Te(t) = 0.7;   % Moderate thermal stress
            Ti(t) = 0.7;   % Moderate training intensity
            Tc(t) = 0.7;   % Moderate complexity
            Tp(t) = 0.7;   % Moderate pressure
            Hr(t) = 0.5;   % Moderate heart rate strain
            Fi(t) = 0.5;   % Moderate fluid intake
            Fl(t) = 0.5;   % Moderate fitness level
            Sh(t) = 0.5;   % Moderate sleep hours
            St(t) = 0.5;   % Moderate sleep timing
        end
    

    elseif case_idx == 4
        % Case 4: Lp=Sp, Lx=0, Fg=0 (Pressure equilibrium, no exhaustion)
        % Good conditions to test if Lp=Sp equilibrium holds without exhaustion
        fprintf('Input Scenario: Good conditions (testing Lp=Sp equilibrium)\n');
        for t = 1:numSteps
            Te(t) = 0;     % No thermal stress
            Ti(t) = 0;     % No training intensity
            Tc(t) = 0;     % No complexity
            Tp(t) = 0;     % No pressure
            Hr(t) = 0;     % No heart rate strain
            Fi(t) = 1;     % Maximum fluid intake
            Fl(t) = 1;     % Maximum fitness level
            Sh(t) = 1;     % Maximum sleep hours
            St(t) = 1;     % Optimal sleep timing
        end

    elseif case_idx == 5
        % Case 5: Lp=1, Lx=1, Fg=Lx (Maximum saturation with Fg tracking Lx)
        % Worst conditions to test maximum saturation equilibrium
        fprintf('Input Scenario: Worst conditions (testing Lp=1, Lx=1, Fg=Lx)\n');
        for t = 1:numSteps
            Te(t) = 1.7;     % Maximum thermal stress
            Ti(t) = 1.7;     % Maximum training intensity
            Tc(t) = 1.7;     % Maximum complexity
            Tp(t) = 1.7;     % Maximum pressure
            Hr(t) = 1.7;     % Maximum heart rate strain
            Fi(t) = 0;     % No fluid intake
            Fl(t) = 0;     % No fitness level
            Sh(t) = 0;     % No sleep hours
            St(t) = 0;     % Poor sleep timing
        end
    end   

    % Initialize state arrays
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
    Lp = zeros(1,numSteps);     % Long-term Experienced Pressure (FORCED)
    Lx = zeros(1,numSteps);     % Long-term Exhaustion (FORCED)
    Fg = zeros(1,numSteps);     % Fatigue (FORCED)
    Rt = zeros(1,numSteps);     % Recovery Time

    % Force initial values
    Lp(1) = Lp_fixed;
    Lx(1) = Lx_fixed;
    Fg(1) = Fg_fixed;
    Sx(1) = Fg_fixed;  % Initialize Sx based on Fg

    %% -------- STEP 1 (t=1) INSTANTANEOUS COMPUTE --------
    % Task Demands
    Td(1) = beta_Td * Ti(1) * (w_Td1*Tc(1) + w_Td2*Tp(1)) + (1 - beta_Td) * (Te(1) * Ti(1));

    % Sleep Quality
    Sq(1) = delta_Sq * Sh(1) + (1-delta_Sq) * St(1);

    % Hydration Level
    Hl(1) = (beta_Hl * Fi(1)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(1) + w_Hl2 * Td(1)));

    % Physical Preparedness
    Pp(1) = beta_Pp * (w_Pp1*Hl(1) + w_Pp2*Fl(1) + w_Pp3*Sq(1)) + (1-beta_Pp) * (1 - Hr(1));

    % Critical Point
    Cp(1) = Pp(1) * (1 - Sx(1));

    % Effort Motivation (depends on Lp - FORCED)
    Em(1) = exp( -eta_Em * ((Td(1) - (1-Lp(1)) ) / alpha_Em) ^ 2 );

    % Generated Effort
    Ge(1) = Em(1) * (w_Ge1*Td(1) + w_Ge2*Cp(1));

    % Recovery Effort
    Re(1) = w_Re1 * min(1, max(0, Cp(1) - Ge(1))) + w_Re2 * Ge(1) + w_Re3 * min(1, max(0, (Pp(1) - Cp(1)) / max(Pp(1), 1e-9)));

    % Short-term Exhaustion (depends on Fg - FORCED)
    Sx(1) = (w_Sx1*Fg(1) + w_Sx2*max(0, Ge(1)-Cp(1))) * (1 - gamma_Sx*Re(1));

    % Short-term Experienced Pressure
    Sp(1) = beta_Sp * (w_Sp1 * Sx(1) + w_Sp2 * max(0, Ge(1)-Cp(1))) + (1-beta_Sp) * (1-Pp(1));

    % Recovery Time
    Rt(1) = (w_Rt1*Sx(1) + w_Rt2*Fg(1)) * (1 - (beta_Rt1*Pp(1) + beta_Rt2*Sq(1) + beta_Rt3*Re(1)));

    %% -------- MAIN LOOP --------
    for t = 2:numSteps
        % Calculate instantaneous variables
        Td(t) = beta_Td * Ti(t) * (w_Td1*Tc(t) + w_Td2*Tp(t)) + (1 - beta_Td) * (Te(t) * Ti(t));
        Sq(t) = delta_Sq * Sh(t) + (1-delta_Sq) * St(t-1);
        Hl(t) = (beta_Hl * Fi(t)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(t) + w_Hl2 * Td(t)));
        Pp(t) = beta_Pp * (w_Pp1*Hl(t) + w_Pp2*Fl(t) + w_Pp3*Sq(t)) + (1-beta_Pp) * (1 - Hr(t));
        Cp(t) = Pp(t) * (1 - Sx(t-1));

        % FORCE temporal variables to equilibrium values
        Lp(t) = Lp_fixed;
        Lx(t) = Lx_fixed;
        Fg(t) = Fg_fixed;

        % Calculate variables that depend on forced values
        Em(t) = exp( -eta_Em * ((Td(t) - (1-Lp(t)) ) / alpha_Em) ^ 2 );
        Ge(t) = Em(t) * (w_Ge1*Td(t) + w_Ge2*Cp(t));
        Re(t) = w_Re1 * min(1, max(0, Cp(t) - Ge(t))) + w_Re2 * Ge(t) + w_Re3 * min(1, max(0, (Pp(t) - Cp(t)) / max(Pp(t), 1e-9)));
        Sx(t) = (w_Sx1*Fg(t) + w_Sx2*max(0, Ge(t)-Cp(t))) * (1 - gamma_Sx*Re(t));
        Sp(t) = beta_Sp * (w_Sp1 * Sx(t-1) + w_Sp2 * max(0, Ge(t)-Cp(t))) + (1-beta_Sp) * (1-Pp(t));
        Rt(t) = (w_Rt1*Sx(t-1) + w_Rt2*Fg(t)) * (1 - (beta_Rt1*Pp(t) + beta_Rt2*Sq(t) + beta_Rt3*Re(t)));
    end

    %% -------- CALCULATE EQUILIBRIUM CONSISTENCY --------
    % Check if computed Sp, Sx are consistent with forced Lp, Lx, Fg

    % Get final values (last 10 time steps average for stability)
    final_idx = numSteps-9:numSteps;

    Sp_final = mean(Sp(final_idx));
    Sx_final = mean(Sx(final_idx));
    Em_final = mean(Em(final_idx));
    Ge_final = mean(Ge(final_idx));
    Cp_final = mean(Cp(final_idx));
    Sx_final = mean(Sx(final_idx));
    Rt_final = mean(Rt(final_idx));

    % Display results
    fprintf('\n--- Final State Values (averaged over last 10 steps) ---\n');
    fprintf('Forced: Lp = %.3f, Lx = %.3f, Fg = %.3f\n', Lp_fixed, Lx_fixed, Fg_fixed);
    fprintf('Computed: Sp = %.3f, Sx = %.3f\n', Sp_final, Sx_final);
    fprintf('          Em = %.3f, Ge = %.3f, Cp = %.3f\n', Em_final, Ge_final, Cp_final);
    fprintf('          Rt = %.3f\n', Rt_final);

    % Check consistency
    fprintf('\n--- Equilibrium Consistency Checks ---\n');

    if case_idx == 3
        % Special checks for Case 3: Lp=Sp, Lx=Sx, Fg=Lx (Moderate Equilibrium)
        fprintf('Testing equilibrium condition: Lp=Sp ∧ Lx=Sx ∧ Fg=Lx\n');
        fprintf('\n');

        % Check 1: Lp = Sp?
        error_Lp_Sp = abs(Lp_fixed - Sp_final);
        fprintf('Check 1 - Lp=Sp: Lp=%.3f, Sp=%.3f, Error=%.3f', Lp_fixed, Sp_final, error_Lp_Sp);
        if error_Lp_Sp < 0.1
            fprintf(' ✓ PASS\n');
        else
            fprintf(' ✗ FAIL (error > 0.1)\n');
        end

        % Check 2: Lx = Sx?
        error_Lx_Sx = abs(Lx_fixed - Sx_final);
        fprintf('Check 2 - Lx=Sx: Lx=%.3f, Sx=%.3f, Error=%.3f', Lx_fixed, Sx_final, error_Lx_Sx);
        if error_Lx_Sx < 0.1
            fprintf(' ✓ PASS\n');
        else
            fprintf(' ✗ FAIL (error > 0.1)\n');
        end

        % Check 3: Fg = Lx?
        error_Fg_Lx = abs(Fg_fixed - Lx_fixed);
        fprintf('Check 3 - Fg=Lx: Fg=%.3f, Lx=%.3f, Error=%.3f', Fg_fixed, Lx_fixed, error_Fg_Lx);
        if error_Fg_Lx < 0.01
            fprintf(' ✓ PASS (forced values match)\n');
        else
            fprintf(' ✗ FAIL\n');
        end

        fprintf('\nInterpretation:\n');
        fprintf('- If Lp≈Sp: Long-term pressure equilibrium is consistent with short-term pressure\n');
        fprintf('- If Lx≈Sx: Long-term exhaustion equilibrium is consistent with short-term exhaustion\n');
        fprintf('- If Fg≈Lx: Fatigue correctly tracks long-term exhaustion\n');

    elseif case_idx == 4
        % Special checks for Case 4: Lp=Sp, Lx=0, Fg=0 (Pressure Equilibrium, No Exhaustion)
        fprintf('Testing equilibrium condition: Lp=Sp (with Lx=0, Fg=0)\n');
        fprintf('\n');

        % Check: Lp = Sp?
        error_Lp_Sp = abs(Lp_fixed - Sp_final);
        fprintf('Check - Lp=Sp: Lp=%.3f, Sp=%.3f, Error=%.3f', Lp_fixed, Sp_final, error_Lp_Sp);
        if error_Lp_Sp < 0.1
            fprintf(' ✓ PASS\n');
        else
            fprintf(' ✗ FAIL (error > 0.1)\n');
        end

        fprintf('\nInterpretation:\n');
        fprintf('- Tests if long-term pressure (Lp=0) equilibrium matches computed short-term pressure\n');
        fprintf('- With no exhaustion (Lx=0, Fg=0) and good conditions, Sp should stabilize at zero\n');

    elseif case_idx == 5
        % Special checks for Case 5: Lp=1, Lx=1, Fg=Lx (Maximum Saturation, Fg tracks Lx)
        fprintf('Testing equilibrium condition: Lp=1 ∧ Lx=1 ∧ Fg=Lx\n');
        fprintf('\n');

        % Check 1: Fg = Lx?
        error_Fg_Lx = abs(Fg_fixed - Lx_fixed);
        fprintf('Check 1 - Fg=Lx: Fg=%.3f, Lx=%.3f, Error=%.3f', Fg_fixed, Lx_fixed, error_Fg_Lx);
        if error_Fg_Lx < 0.01
            fprintf(' ✓ PASS (forced values match)\n');
        else
            fprintf(' ✗ FAIL\n');
        end

        % Check 2: Lx = Sx?
        error_Lx_Sx = abs(Lx_fixed - Sx_final);
        fprintf('Check 2 - Lx=Sx: Lx=%.3f, Sx=%.3f, Error=%.3f', Lx_fixed, Sx_final, error_Lx_Sx);
        if error_Lx_Sx < 0.1
            fprintf(' ✓ PASS\n');
        else
            fprintf(' ✗ FAIL (error > 0.1)\n');
        end

        fprintf('\nInterpretation:\n');
        fprintf('- Tests maximum saturation: all temporal variables at 1.0\n');
        fprintf('- Fg=Lx verifies fatigue correctly tracks exhaustion at maximum\n');
        fprintf('- Lx≈Sx verifies long-term exhaustion equilibrium with short-term exhaustion\n');

    else
        % Standard checks for other cases
        % Check 1: If Lp is forced, does computed Sp make sense?
        if Lp_fixed == 0
            fprintf('Lp=0 forced → Expected Sp ≈ 0 (if no pressure inputs)\n');
            fprintf('             Actual Sp = %.3f\n', Sp_final);
        elseif Lp_fixed == 1
            fprintf('Lp=1 forced → Expected Sp ≈ 1 (saturated pressure)\n');
            fprintf('             Actual Sp = %.3f\n', Sp_final);
        end

        % Check 2: If Fg is forced, does Sx make sense?
        fprintf('Fg=%.1f forced → Sx influenced by w_Sx1*Fg = %.3f\n', Fg_fixed, w_Sx1*Fg_fixed);
        fprintf('             Actual Sx = %.3f\n', Sx_final);

        % Check 3: If Lp affects Em
        expected_Em_arg = (Td(end) - (1-Lp_fixed)) / alpha_Em;
        fprintf('Lp=%.1f forced → Em argument = (%.3f - %.3f)/%.2f = %.3f\n', ...
                Lp_fixed, Td(end), (1-Lp_fixed), alpha_Em, expected_Em_arg);
        fprintf('             Actual Em = %.3f\n', Em_final);
    end

    % Store results for plotting
    results{case_idx}.name = test_cases(case_idx).name;
    results{case_idx}.time = time;
    results{case_idx}.Lp = Lp;
    results{case_idx}.Lx = Lx;
    results{case_idx}.Fg = Fg;
    results{case_idx}.Sp = Sp;
    results{case_idx}.Sx = Sx;
    results{case_idx}.Em = Em;
    results{case_idx}.Ge = Ge;
    results{case_idx}.Cp = Cp;
    results{case_idx}.Re = Re;
    results{case_idx}.Rt = Rt;
    results{case_idx}.Pp = Pp;

    %% -------- CREATE SEPARATE PLOTS FOR THIS CASE --------
    minY = -0.05;
    maxY = 1.05;

    % Figure 1: Forced Temporal Variables (Lp, Lx, Fg)
    fig_name = sprintf('%s - Forced Temporal Variables', test_cases(case_idx).name);
    figure('Name', fig_name, 'Position', [100, 100, 1200, 800]);

    subplot(3,1,1);
    plot(time, Lp, 'LineWidth', 2, 'Color', [0.00, 0.45, 0.74]);
    ylabel('Lp (Long-term Pressure)');
    title(sprintf('%s: Forced Temporal Variables', test_cases(case_idx).name));
    ylim([minY maxY]);
    grid on;

    subplot(3,1,2);
    plot(time, Lx, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.10]);
    ylabel('Lx (Long-term Exhaustion)');
    ylim([minY maxY]);
    grid on;

    subplot(3,1,3);
    plot(time, Fg, 'LineWidth', 2, 'Color', [0.93, 0.69, 0.13]);
    xlabel('Time');
    ylabel('Fg (Fatigue)');
    ylim([minY maxY]);
    grid on;

    % Figure 2: Computed Short-term Variables (Sp, Sx)
    fig_name = sprintf('%s - Short-term Variables', test_cases(case_idx).name);
    figure('Name', fig_name, 'Position', [150, 150, 1200, 600]);

    subplot(2,1,1);
    plot(time, Sp, 'LineWidth', 2, 'Color', [0.49, 0.18, 0.56]);
    ylabel('Sp (Short-term Pressure)');
    title(sprintf('%s: Computed Short-term Variables', test_cases(case_idx).name));
    ylim([minY maxY]);
    grid on;

    subplot(2,1,2);
    plot(time, Sx, 'LineWidth', 2, 'Color', [0.47, 0.67, 0.19]);
    xlabel('Time');
    ylabel('Sx (Short-term Exhaustion)');
    ylim([minY maxY]);
    grid on;

    % Figure 3: Recovery Time and Influencing Variables (Sx, Fg, Rt)
    fig_name = sprintf('%s - Recovery Time Analysis', test_cases(case_idx).name);
    figure('Name', fig_name, 'Position', [200, 200, 1200, 800]);

    subplot(3,1,1);
    plot(time, Sx, 'LineWidth', 2, 'Color', [0.47, 0.67, 0.19]);
    ylabel('Sx (Short-term Exhaustion)');
    title(sprintf('%s: Recovery Time and Influencing Variables', test_cases(case_idx).name));
    ylim([minY maxY]);
    grid on;

    subplot(3,1,2);
    plot(time, Fg, 'LineWidth', 2, 'Color', [0.93, 0.69, 0.13]);
    ylabel('Fg (Fatigue)');
    ylim([minY maxY]);
    grid on;

    subplot(3,1,3);
    plot(time, Rt, 'LineWidth', 2, 'Color', [0.00, 0.60, 0.50]);
    xlabel('Time');
    ylabel('Rt (Recovery Time)');
    ylim([minY maxY]);
    grid on;

    % Figure 4: Effort Dynamics (Em, Ge, Cp)
    fig_name = sprintf('%s - Effort Dynamics', test_cases(case_idx).name);
    figure('Name', fig_name, 'Position', [250, 250, 1200, 800]);

    subplot(3,1,1);
    plot(time, Em, 'LineWidth', 2, 'Color', [0.64, 0.08, 0.18]);
    ylabel('Em (Effort Motivation)');
    title(sprintf('%s: Effort Dynamics', test_cases(case_idx).name));
    ylim([minY maxY]);
    grid on;

    subplot(3,1,2);
    plot(time, Ge, 'LineWidth', 2, 'Color', [0.30, 0.75, 0.93]);
    ylabel('Ge (Generated Effort)');
    ylim([minY maxY]);
    grid on;

    subplot(3,1,3);
    plot(time, Cp, 'LineWidth', 2, 'Color', [0.47, 0.25, 0.80]);
    xlabel('Time');
    ylabel('Cp (Critical Point)');
    ylim([minY maxY]);
    grid on;
end

fprintf('\n========================================\n');
fprintf('All test cases completed!\n');
fprintf('Total figures created: %d (4 per case)\n', num_cases * 4);
fprintf('========================================\n');
