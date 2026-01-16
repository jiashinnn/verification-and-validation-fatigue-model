% ======================================================================
%  VALIDATION STEP 2: SENSITIVITY ANALYSIS - INPUT VARIABLES (OAT)
%  Octave/MATLAB compatible (single-file, runs from Command Window)
%
%  HOW TO RUN:
%   1) Save this file as:  validation_step2_sensitivity_analysis_inputs.m
%   2) In Octave/MATLAB Command Window (same folder):
%         validation_step2_sensitivity_analysis_inputs
% ======================================================================

function validation_step2_sensitivity_analysis_inputs()

clc; clear; close all;

fprintf('=== SENSITIVITY ANALYSIS: INPUT VARIABLES ===\n\n');

%% ===================== SIMULATION PARAMETERS =====================
T  = 50;           % total duration
dt = 0.1;          % step size
time = 0:dt:T;     % exact dt steps
numSteps = numel(time);

%% ===================== MODEL PARAMETERS (FROM YOUR MODEL) =====================
% Task Demands
w_Td1 = 0.4; w_Td2 = 0.6; beta_Td = 0.7;

% Sleep Quality
delta_Sq = 0.7;

% Hydration Level
beta_Hl = 0.5; w_Hl1 = 0.6; w_Hl2 = 0.4;

% Physical Preparedness
beta_Pp  = 0.7; w_Pp1 = 0.25; w_Pp2 = 0.30; w_Pp3 = 0.45;

% Effort Motivation
alpha_Em = 0.4; eta_Em = 0.5;

% Generated Effort
w_Ge1 = 0.85; w_Ge2 = 0.15;

% Recovery Effort (CALIBRATED)
w_Re1 = 0.25; w_Re2 = 0.20; w_Re3 = 0.20;

% Short-term Exhaustion (CALIBRATED)
w_Sx1 = 0.70; w_Sx2 = 1.10; gamma_Sx = 0.50;

% Short-term Experienced Pressure
w_Sp1 = 0.6; w_Sp2 = 0.4; beta_Sp = 0.7;

% Recovery Time
w_Rt1 = 0.5; w_Rt2 = 0.5;
beta_Rt1 = 0.3; beta_Rt2 = 0.3; beta_Rt3 = 0.4;

%% ===================== MULTI-POINT BASELINE VALUES =====================
% Two baseline operating points to test sensitivity at different regions
baseline_values = [0.5, 0.8];
baseline_labels = {'Nominal (0.5)', 'High (0.8)'};
n_baselines = length(baseline_values);

%% ===================== SENSITIVITY TEST VALUES =====================
test_values = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
n_tests = length(test_values);

% Input names
input_names = {'Te', 'Ti', 'Tc', 'Tp', 'Fi', 'Fl', 'Hr', 'Sh', 'St'};
input_labels = {'Thermal Env', 'Training Intensity', 'Training Complexity', ...
                'Training Pressure', 'Fluid Intake', 'Fitness Level', ...
                'Heart Rate', 'Sleep Hours', 'Sleep Timing'};
n_inputs = length(input_names);

% Store results for each baseline
sensitivity_results = struct();
for b = 1:n_baselines
    baseline_name = sprintf('baseline_%d', round(baseline_values(b)*10));
    sensitivity_results.(baseline_name) = struct();
    for i = 1:n_inputs
        sensitivity_results.(baseline_name).(input_names{i}) = zeros(n_tests, numSteps);
    end
end

%% ===================== RUN MULTI-POINT OAT SENSITIVITY ANALYSIS =====================
fprintf('Running Multi-Point One-At-a-Time (OAT) Sensitivity Analysis...\n');
fprintf('Testing %d baselines × %d inputs × %d values = %d simulations\n\n', ...
        n_baselines, n_inputs, n_tests, n_baselines*n_inputs*n_tests);

for baseline_idx = 1:n_baselines
    baseline_val = baseline_values(baseline_idx);
    baseline_name = sprintf('baseline_%d', round(baseline_val*10));

    fprintf('=== BASELINE %d/%d: All other inputs = %.1f ===\n', ...
            baseline_idx, n_baselines, baseline_val);

    for input_idx = 1:n_inputs
        current_input = input_names{input_idx};
        fprintf('  [%d/%d] Testing %s...\n', input_idx, n_inputs, current_input);

        for test_idx = 1:n_tests
            test_val = test_values(test_idx);

            % Set all inputs to current baseline (constant over time)
            Te = ones(1, numSteps) * baseline_val;
            Ti = ones(1, numSteps) * baseline_val;
            Tc = ones(1, numSteps) * baseline_val;
            Tp = ones(1, numSteps) * baseline_val;
            Fi = ones(1, numSteps) * baseline_val;
            Fl = ones(1, numSteps) * baseline_val;
            Hr = ones(1, numSteps) * baseline_val;
            Sh = ones(1, numSteps) * baseline_val;
            St = ones(1, numSteps) * baseline_val;

            % Override the current input being tested
            switch current_input
                case 'Te', Te = ones(1, numSteps) * test_val;
                case 'Ti', Ti = ones(1, numSteps) * test_val;
                case 'Tc', Tc = ones(1, numSteps) * test_val;
                case 'Tp', Tp = ones(1, numSteps) * test_val;
                case 'Fi', Fi = ones(1, numSteps) * test_val;
                case 'Fl', Fl = ones(1, numSteps) * test_val;
                case 'Hr', Hr = ones(1, numSteps) * test_val;
                case 'Sh', Sh = ones(1, numSteps) * test_val;
                case 'St', St = ones(1, numSteps) * test_val;
            end

            % Run simulation
            Fg = run_fatigue_simulation(Te, Ti, Tc, Tp, Fi, Fl, Hr, Sh, St, ...
                                         numSteps, dt, ...
                                         w_Td1, w_Td2, beta_Td, delta_Sq, beta_Hl, w_Hl1, w_Hl2, ...
                                         beta_Pp, w_Pp1, w_Pp2, w_Pp3, alpha_Em, eta_Em, ...
                                         w_Ge1, w_Ge2, w_Re1, w_Re2, w_Re3, ...
                                         w_Sx1, w_Sx2, gamma_Sx, w_Sp1, w_Sp2, beta_Sp, ...
                                         w_Rt1, w_Rt2, beta_Rt1, beta_Rt2, beta_Rt3);

            sensitivity_results.(baseline_name).(current_input)(test_idx, :) = Fg;
        end
    end
    fprintf('\n');
end

fprintf('✓ Multi-point sensitivity analysis complete!\n\n');

%% ===================== VISUALIZATION =====================
fprintf('Creating visualizations...\n');

colors = lines(n_tests);

% Create one figure for each baseline
for baseline_idx = 1:n_baselines
    baseline_val = baseline_values(baseline_idx);
    baseline_name = sprintf('baseline_%d', round(baseline_val*10));

    figure('Name', sprintf('Sensitivity Analysis - Baseline %.1f', baseline_val), ...
           'Position', [100 + (baseline_idx-1)*50, 100 + (baseline_idx-1)*50, 1400, 900]);

    for i = 1:n_inputs
        subplot(3, 3, i);
        hold on;

        for test_idx = 1:n_tests
            plot(time, sensitivity_results.(baseline_name).(input_names{i})(test_idx, :), ...
                 'LineWidth', 2, 'Color', colors(test_idx, :), ...
                 'DisplayName', sprintf('%s=%.1f', input_names{i}, test_values(test_idx)));
        end

        xlabel('Time (t)');
        ylabel('Fatigue (Fg)');
        title(sprintf('Sensitivity to %s', input_labels{i}));
        grid on;
        ylim([0, 1]);

        legend('Location', 'northeast');
        set(gca, 'FontSize', 9);
    end

    % ---- Octave-safe overall title (sgtitle replacement) ----
    title_str = sprintf('OAT Input Sensitivity Analysis | Baseline: All Other Inputs = %.1f', baseline_val);
    if exist('sgtitle','file') == 2
        sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    else
        % Create a figure-level title using an annotation textbox
        annotation(gcf,'textbox',[0 0.96 1 0.04], 'String', title_str, ...
            'EdgeColor','none', 'HorizontalAlignment','center', ...
            'FontSize', 14, 'FontWeight','bold');
    end
end


%% ===================== CALCULATE SENSITIVITY INDICES =====================
fprintf('--- SENSITIVITY INDICES (Final Fatigue at t=%.1f) ---\n', T);

sensitivity_indices = struct();

for baseline_idx = 1:n_baselines
    baseline_val = baseline_values(baseline_idx);
    baseline_name = sprintf('baseline_%d', round(baseline_val*10));

    fprintf('\n>>> Baseline %.1f (All other inputs = %.1f):\n', baseline_val, baseline_val);

    sensitivity_indices.(baseline_name) = struct();
    ranges = zeros(1, n_inputs);

    for i = 1:n_inputs
        final_fatigues = sensitivity_results.(baseline_name).(input_names{i})(:, end);

        sensitivity_indices.(baseline_name).(input_names{i}) = struct();
        sensitivity_indices.(baseline_name).(input_names{i}).values = test_values';
        sensitivity_indices.(baseline_name).(input_names{i}).final_fatigue = final_fatigues;
        sensitivity_indices.(baseline_name).(input_names{i}).range = max(final_fatigues) - min(final_fatigues);
        sensitivity_indices.(baseline_name).(input_names{i}).mean = mean(final_fatigues);

        ranges(i) = sensitivity_indices.(baseline_name).(input_names{i}).range;

        fprintf('  %-20s Range: %.4f  Mean: %.4f\n', input_labels{i}, ...
                sensitivity_indices.(baseline_name).(input_names{i}).range, ...
                sensitivity_indices.(baseline_name).(input_names{i}).mean);
    end

    % Store ranges for bar chart
    sensitivity_indices.(baseline_name).all_ranges = ranges;
end

%% ===================== BAR CHART OF RANGES (All Baselines) =====================
fprintf('\nCreating comparison bar chart...\n');

figure('Name', 'Multi-Point Sensitivity Comparison', 'Position', [200, 200, 1400, 600]);

% Prepare data matrix: rows = baselines, columns = inputs
range_matrix = zeros(n_baselines, n_inputs);
for baseline_idx = 1:n_baselines
    baseline_name = sprintf('baseline_%d', round(baseline_values(baseline_idx)*10));
    range_matrix(baseline_idx, :) = sensitivity_indices.(baseline_name).all_ranges;
end

% Create grouped bar chart
bar_handle = bar(range_matrix');
set(gca, 'XTick', 1:n_inputs, 'XTickLabel', input_labels, 'XTickLabelRotation', 35);
ylabel('Fatigue Range (Max - Min) at Final Time');
title('Multi-Point OAT Input Sensitivity Comparison');
legend(baseline_labels, 'Location', 'northeast');
grid on;

%% ===================== SAVE RESULTS =====================
save('sensitivity_analysis_inputs_multipoint.mat', 'sensitivity_results', 'sensitivity_indices', ...
     'test_values', 'input_names', 'input_labels', 'time', 'T', 'dt', ...
     'baseline_values', 'baseline_labels', 'n_baselines');

fprintf('\n✓ Results saved: sensitivity_analysis_inputs_multipoint.mat\n');
fprintf('✓ Visualizations created: %d figures (%d baselines + 1 comparison)\n\n', n_baselines + 1, n_baselines);

end % end main function


%% ===================== HELPER FUNCTION =====================
function Fg = run_fatigue_simulation(Te, Ti, Tc, Tp, Fi, Fl, Hr, Sh, St, ...
                                      numSteps, dt, ...
                                      w_Td1, w_Td2, beta_Td, delta_Sq, beta_Hl, w_Hl1, w_Hl2, ...
                                      beta_Pp, w_Pp1, w_Pp2, w_Pp3, alpha_Em, eta_Em, ...
                                      w_Ge1, w_Ge2, w_Re1, w_Re2, w_Re3, ...
                                      w_Sx1, w_Sx2, gamma_Sx, w_Sp1, w_Sp2, beta_Sp, ...
                                      w_Rt1, w_Rt2, beta_Rt1, beta_Rt2, beta_Rt3)

    % Initialize state arrays
    Td = zeros(1, numSteps); Hl = zeros(1, numSteps); Sq = zeros(1, numSteps);
    Pp = zeros(1, numSteps); Em = zeros(1, numSteps); Ge = zeros(1, numSteps);
    Cp = zeros(1, numSteps); Sp = zeros(1, numSteps); Sx = zeros(1, numSteps);
    Re = zeros(1, numSteps); Lp = zeros(1, numSteps); Lx = zeros(1, numSteps);
    Fg = zeros(1, numSteps); Rt = zeros(1, numSteps);

    % Initial values
    Lp(1) = 0.20; Lx(1) = 0.20; Fg(1) = 0.20; Sx(1) = 0.20;

    % Step 1
    Td(1) = beta_Td * Ti(1) * (w_Td1*Tc(1) + w_Td2*Tp(1)) + (1 - beta_Td) * (Te(1) * Ti(1));
    Sq(1) = delta_Sq * Sh(1) + (1-delta_Sq) * St(1);
    Hl(1) = (beta_Hl * Fi(1)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(1) + w_Hl2 * Td(1)));
    Pp(1) = beta_Pp * (w_Pp1*Hl(1) + w_Pp2*Fl(1) + w_Pp3*Sq(1)) + (1-beta_Pp) * (1 - Hr(1));
    Cp(1) = Pp(1) * (1 - Sx(1));
    Em(1) = exp( -eta_Em * ((Td(1) - (1-Lp(1)) ) / alpha_Em) ^ 2 );
    Ge(1) = Em(1) * (w_Ge1*Td(1) + w_Ge2*Cp(1));
    Re(1) = w_Re1 * min(1, max(0, Cp(1) - Ge(1))) ...
          + w_Re2 * Ge(1) ...
          + w_Re3 * min(1, max(0, (Pp(1) - Cp(1)) / max(Pp(1), 1e-9)));
    Sx(1) = (w_Sx1*Fg(1) + w_Sx2*max(0, Ge(1)-Cp(1))) * (1 - gamma_Sx*Re(1));
    Sp(1) = beta_Sp * (w_Sp1 * Sx(1) + w_Sp2 * max(0, Ge(1)-Cp(1))) + (1-beta_Sp) * (1-Pp(1));
    Rt(1) = (w_Rt1*Sx(1) + w_Rt2*Fg(1)) * (1 - (beta_Rt1*Pp(1) + beta_Rt2*Sq(1) + beta_Rt3*Re(1)));

    % Main loop
    for t = 2:numSteps
        Td(t) = beta_Td * Ti(t) * (w_Td1*Tc(t) + w_Td2*Tp(t)) + (1 - beta_Td) * (Te(t) * Ti(t));
        Sq(t) = delta_Sq * Sh(t) + (1-delta_Sq) * St(t-1);
        Hl(t) = (beta_Hl * Fi(t)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(t) + w_Hl2 * Td(t)));
        Pp(t) = beta_Pp * (w_Pp1*Hl(t) + w_Pp2*Fl(t) + w_Pp3*Sq(t)) + (1-beta_Pp) * (1 - Hr(t));
        Cp(t) = Pp(t) * (1 - Sx(t-1));
        Em(t) = exp( -eta_Em * ((Td(t) - (1-Lp(t-1)) ) / alpha_Em) ^ 2 );
        Ge(t) = Em(t) * (w_Ge1*Td(t) + w_Ge2*Cp(t));
        Re(t) = w_Re1 * min(1, max(0, Cp(t) - Ge(t))) ...
              + w_Re2 * Ge(t) ...
              + w_Re3 * min(1, max(0, (Pp(t) - Cp(t)) / max(Pp(t), 1e-9)));
        Sx(t) = (w_Sx1*Fg(t-1) + w_Sx2*max(0, Ge(t)-Cp(t))) * (1 - gamma_Sx*Re(t));
        Sp(t) = beta_Sp * (w_Sp1 * Sx(t-1) + w_Sp2 * max(0, Ge(t)-Cp(t))) + (1-beta_Sp) * (1-Pp(t));

        Lp(t) = Lp(t-1) + (Sp(t) - Lp(t-1)) * Lp(t-1) * (1-Lp(t-1)) * dt;
        Lx(t) = Lx(t-1) + (Sx(t-1) - Lx(t-1)) * Lx(t-1) * (1-Lx(t-1)) * dt;
        Fg(t) = Fg(t-1) + (Lx(t-1) - Fg(t-1)) * Fg(t-1) * (1-Fg(t-1)) * dt;

        Rt(t) = (w_Rt1*Sx(t-1) + w_Rt2*Fg(t-1)) * (1 - (beta_Rt1*Pp(t) + beta_Rt2*Sq(t) + beta_Rt3*Re(t)));
    end

end

