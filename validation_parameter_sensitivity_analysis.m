% ======================================================================
%  VALIDATION STEP 3: COMPREHENSIVE PARAMETER SENSITIVITY ANALYSIS
%  Purpose: Test sensitivity to ALL model parameters
%  Method: One-at-a-time testing with extreme input scenario
%  Display Order: Priority based on proximity to Fg equation
%
%  HOW TO RUN:
%   1) Save this file as: validation_step3_sensitivity_parameters_complete.m
%   2) In MATLAB Command Window (same folder):
%         validation_step3_sensitivity_parameters_complete
% ======================================================================

function validation_step3_sensitivity_parameters_complete()

clc; clear; close all;

fprintf('=== COMPREHENSIVE PARAMETER SENSITIVITY ANALYSIS ===\n\n');

%% ===================== SIMULATION PARAMETERS =====================
T  = 50;
dt = 0.1;
time = 0:dt:T;
numSteps = numel(time);

%% ===================== EXTREME INPUT SCENARIO =====================
% Critical worst-case scenario
fprintf('Using EXTREME INPUT SCENARIO (worst case):\n');
fprintf('  High stress:  Te=Ti=Tc=Tp=Hr = 1.0\n');
fprintf('  Poor recovery: Sh=St=Fi=Fl = 0.0\n\n');

Te_fixed = ones(1, numSteps) * 1.0;
Ti_fixed = ones(1, numSteps) * 0.8;
Tc_fixed = ones(1, numSteps) * 0.8;
Tp_fixed = ones(1, numSteps) * 0.8;
Hr_fixed = ones(1, numSteps) * 0.75;
Sh_fixed = ones(1, numSteps) * 0.55;
St_fixed = ones(1, numSteps) * 0.4;
Fi_fixed = ones(1, numSteps) * 0.4;
Fl_fixed = ones(1, numSteps) * 0.4;

%% ===================== BASELINE PARAMETER VALUES (CALIBRATED) =====================
% From refined_model.m
params_baseline = struct();
params_baseline.w_Td1 = 0.4; params_baseline.w_Td2 = 0.6; params_baseline.beta_Td = 0.7;
params_baseline.delta_Sq = 0.7;
params_baseline.beta_Hl = 0.5; params_baseline.w_Hl1 = 0.6; params_baseline.w_Hl2 = 0.4;
params_baseline.beta_Pp = 0.7; params_baseline.w_Pp1 = 0.25; params_baseline.w_Pp2 = 0.30; params_baseline.w_Pp3 = 0.45;
params_baseline.alpha_Em = 0.4; params_baseline.eta_Em = 0.5;
params_baseline.w_Ge1 = 0.85; params_baseline.w_Ge2 = 0.15;
params_baseline.w_Re1 = 0.25; params_baseline.w_Re2 = 0.20; params_baseline.w_Re3 = 0.20;  % CALIBRATED
params_baseline.w_Sx1 = 0.70; params_baseline.w_Sx2 = 1.10; params_baseline.gamma_Sx = 0.50;  % CALIBRATED
params_baseline.w_Sp1 = 0.6; params_baseline.w_Sp2 = 0.4; params_baseline.beta_Sp = 0.7;
params_baseline.w_Rt1 = 0.5; params_baseline.w_Rt2 = 0.5;
params_baseline.beta_Rt1 = 0.3; params_baseline.beta_Rt2 = 0.3; params_baseline.beta_Rt3 = 0.4;

%% ===================== DEFINE PARAMETER TESTS (IN PRIORITY ORDER) =====================
% Test values
test_independent = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
test_independent_alpha = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0];  % Special for alpha_Em (finer resolution at low values)
test_weights = [0.1, 0.3, 0.5, 0.7, 0.9];

% Parameter definitions (18 items in priority order)
param_tests = {};

% PRIORITY 1: Short-term Exhaustion (closest to Fg)
param_tests{end+1} = struct('name', 'w_Sx1', 'label', 'w_{Sx1} (Short exhaustion)', ...
    'type', 'independent', 'values', test_independent, 'priority', 1);
param_tests{end+1} = struct('name', 'w_Sx2', 'label', 'w_{Sx2} (Short exhaustion)', ...
    'type', 'independent', 'values', test_independent, 'priority', 1);
param_tests{end+1} = struct('name', 'gamma_Sx', 'label', '\gamma_{Sx} (Exhaustion recovery)', ...
    'type', 'independent', 'values', test_independent, 'priority', 1);

% PRIORITY 2: Generated Effort & Recovery Effort
param_tests{end+1} = struct('name', 'w_Ge', 'label', 'w_{Ge1}/w_{Ge2} (Generated effort weights)', ...
    'type', 'pair', 'param1', 'w_Ge1', 'param2', 'w_Ge2', 'values', test_weights, 'priority', 2);
param_tests{end+1} = struct('name', 'w_Re1', 'label', 'w_{Re1} (Recovery effort)', ...
    'type', 'triplet', 'param1', 'w_Re1', 'param2', 'w_Re2', 'param3', 'w_Re3', ...
    'values', test_weights, 'priority', 2);

% PRIORITY 3: Short-term Pressure
param_tests{end+1} = struct('name', 'beta_Sp', 'label', '\beta_{Sp} (Experienced pressure weight)', ...
    'type', 'independent', 'values', test_independent, 'priority', 3);
param_tests{end+1} = struct('name', 'w_Sp', 'label', 'w_{Sp1}/w_{Sp2} (Experienced pressure weights)', ...
    'type', 'pair', 'param1', 'w_Sp1', 'param2', 'w_Sp2', 'values', test_weights, 'priority', 3);

% PRIORITY 4: Physical Preparedness
param_tests{end+1} = struct('name', 'beta_Pp', 'label', '\beta_{Pp} (Physical preparedness weight)', ...
    'type', 'independent', 'values', test_independent, 'priority', 4);
param_tests{end+1} = struct('name', 'w_Pp1', 'label', 'w_{Pp1} (Physical preparedness weights)', ...
    'type', 'triplet', 'param1', 'w_Pp1', 'param2', 'w_Pp2', 'param3', 'w_Pp3', ...
    'values', test_weights, 'priority', 4);

% PRIORITY 5: Effort Motivation
param_tests{end+1} = struct('name', 'alpha_Em', 'label', '\alpha_{Em} (Effort Motivation)', ...
    'type', 'independent', 'values', test_independent_alpha, 'priority', 5);
param_tests{end+1} = struct('name', 'eta_Em', 'label', '\eta_{Em} (Effort Motivation)', ...
    'type', 'independent', 'values', test_independent, 'priority', 5);

% PRIORITY 6: Sleep Quality & Hydration
param_tests{end+1} = struct('name', 'delta_Sq', 'label', '\delta_{Sq} (Sleep quality balance)', ...
    'type', 'independent', 'values', test_independent, 'priority', 6);
param_tests{end+1} = struct('name', 'beta_Hl', 'label', '\beta_{Hl} (Hydration level weight)', ...
    'type', 'independent', 'values', test_independent, 'priority', 6);
param_tests{end+1} = struct('name', 'w_Hl', 'label', 'w_{Hl1}/w_{Hl2} (Hydration level weights)', ...
    'type', 'pair', 'param1', 'w_Hl1', 'param2', 'w_Hl2', 'values', test_weights, 'priority', 6);

% PRIORITY 7: Task Demands
param_tests{end+1} = struct('name', 'w_Td', 'label', 'w_{Td1}/w_{Td2} (Task demand weights)', ...
    'type', 'pair', 'param1', 'w_Td1', 'param2', 'w_Td2', 'values', test_weights, 'priority', 7);

% PRIORITY 8: Recovery Time (lowest priority, doesn't affect Fg)
param_tests{end+1} = struct('name', 'w_Rt', 'label', 'w_{Rt1}/w_{Rt2} (Recovery time weights)', ...
    'type', 'pair', 'param1', 'w_Rt1', 'param2', 'w_Rt2', 'values', test_weights, 'priority', 8);
param_tests{end+1} = struct('name', 'beta_Rt1', 'label', '\beta_{Rt1} (Recovery time factors)', ...
    'type', 'triplet', 'param1', 'beta_Rt1', 'param2', 'beta_Rt2', 'param3', 'beta_Rt3', ...
    'values', test_weights, 'priority', 8);

n_params = length(param_tests);

%% ===================== RUN PARAMETER SENSITIVITY TESTS =====================
fprintf('Running parameter sensitivity tests...\n');
fprintf('Total parameter items to test: %d\n\n', n_params);

param_results = cell(n_params, 1);

for i = 1:n_params
    ptest = param_tests{i};
    fprintf('[%2d/%2d] Testing: %s (Priority %d)...\n', i, n_params, ptest.name, ptest.priority);

    n_values = length(ptest.values);
    results = zeros(n_values, numSteps);

    for j = 1:n_values
        % Start with baseline parameters
        p = params_baseline;

        % Modify parameter(s) based on type
        switch ptest.type
            case 'independent'
                % Single parameter variation
                p.(ptest.name) = ptest.values(j);

            case 'pair'
                % Two parameters that sum to 1
                p.(ptest.param1) = ptest.values(j);
                p.(ptest.param2) = 1.0 - ptest.values(j);

            case 'triplet'
                % Three parameters, vary first, split rest equally
                p.(ptest.param1) = ptest.values(j);
                remaining = (1.0 - ptest.values(j)) / 2.0;
                p.(ptest.param2) = remaining;
                p.(ptest.param3) = remaining;
        end

        % Run simulation with modified parameters
        Fg = run_simulation(Te_fixed, Ti_fixed, Tc_fixed, Tp_fixed, ...
                            Fi_fixed, Fl_fixed, Hr_fixed, Sh_fixed, St_fixed, ...
                            numSteps, dt, p);

        results(j, :) = Fg;
    end

    param_results{i} = results;
end

fprintf('\n✓ Parameter sensitivity complete!\n\n');

%% ===================== CALCULATE SENSITIVITY INDICES =====================
fprintf('--- PARAMETER SENSITIVITY SUMMARY (Priority Order) ---\n');
fprintf('%-40s Priority  Range   Mean\n', 'Parameter');
fprintf('%s\n', repmat('-', 1, 75));

sensitivity_summary = struct();
for i = 1:n_params
    ptest = param_tests{i};
    final_fatigues = param_results{i}(:, end);

    range_val = max(final_fatigues) - min(final_fatigues);
    mean_val = mean(final_fatigues);

    sensitivity_summary(i).name = ptest.name;
    sensitivity_summary(i).label = ptest.label;
    sensitivity_summary(i).priority = ptest.priority;
    sensitivity_summary(i).range = range_val;
    sensitivity_summary(i).mean = mean_val;

    fprintf('%-40s %d        %.4f  %.4f\n', ptest.label, ptest.priority, range_val, mean_val);
end
fprintf('\n');

%% ===================== VISUALIZATION =====================
fprintf('Creating visualizations...\n');

% Determine grid layout for subplots
n_cols = 4;
n_rows = ceil(n_params / n_cols);

% Figure: Parameter Sensitivity in Priority Order
fig1 = figure('Name', 'Parameter Sensitivity Analysis (Priority Order)', ...
              'Position', [50, 50, 1600, 1200]);

colors = lines(6); % Color palette for test values

for i = 1:n_params
    subplot(n_rows, n_cols, i);
    hold on;

    ptest = param_tests{i};
    n_values = length(ptest.values);

    for j = 1:n_values
        plot(time, param_results{i}(j, :), 'LineWidth', 1.5, ...
             'Color', colors(mod(j-1, 6)+1, :), ...
             'DisplayName', sprintf('%.2f', ptest.values(j)));
    end

    xlabel('Time', 'FontSize', 8);
    ylabel('Fatigue (Fg)', 'FontSize', 8);
    title(sprintf('[P%d] %s', ptest.priority, ptest.label), 'FontSize', 9, 'Interpreter', 'tex');
    legend('Location', 'best', 'FontSize', 7);
    grid on;
    ylim([0, 1]);
    xlim([0, T]);
end

% Add overall title
annotation(fig1, 'textbox', [0 0.96 1 0.04], ...
    'String', 'Parameter Sensitivity Analysis (Ordered by Proximity to Fg Equation)', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ===================== BAR CHART: SENSITIVITY RANKING =====================
fig2 = figure('Name', 'Parameter Sensitivity Ranking', 'Position', [100, 100, 1200, 600]);

ranges = [sensitivity_summary.range];
param_labels = {sensitivity_summary.label};
priorities = [sensitivity_summary.priority];

% Color by priority
priority_colors = [0.8 0.2 0.2;   % Priority 1: Red
                   1.0 0.5 0.2;   % Priority 2: Orange
                   1.0 0.8 0.2;   % Priority 3: Yellow
                   0.6 0.8 0.3;   % Priority 4: Yellow-green
                   0.3 0.8 0.6;   % Priority 5: Green
                   0.2 0.6 0.8;   % Priority 6: Blue
                   0.5 0.4 0.8;   % Priority 7: Purple
                   0.7 0.5 0.7];  % Priority 8: Violet

bar_colors = zeros(n_params, 3);
for i = 1:n_params
    bar_colors(i, :) = priority_colors(priorities(i), :);
end

b = barh(1:n_params, ranges);
b.FaceColor = 'flat';
b.CData = bar_colors;

set(gca, 'YTick', 1:n_params, 'YTickLabel', param_labels, 'FontSize', 8);
set(gca, 'TickLabelInterpreter', 'tex');
xlabel('Sensitivity Range (Max - Min Final Fatigue)', 'FontSize', 12);
title('Parameter Sensitivity Ranking (Priority Order: P1 Top → P8 Bottom)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add value labels
for i = 1:n_params
    text(ranges(i) + 0.01, i, sprintf('%.4f', ranges(i)), ...
         'FontSize', 9, 'VerticalAlignment', 'middle');
end

%% ===================== FIGURES BY PRIORITY GROUP =====================
fprintf('Creating priority-grouped figures...\n');

% Group parameters by priority
for priority_level = 1:8
    % Find all parameters with this priority
    priority_indices = find([sensitivity_summary.priority] == priority_level);

    if isempty(priority_indices)
        continue; % Skip if no parameters at this priority
    end

    n_params_priority = length(priority_indices);

    % Determine subplot layout
    if n_params_priority == 1
        n_rows_p = 1; n_cols_p = 1;
    elseif n_params_priority == 2
        n_rows_p = 1; n_cols_p = 2;
    elseif n_params_priority <= 4
        n_rows_p = 2; n_cols_p = 2;
    elseif n_params_priority <= 6
        n_rows_p = 2; n_cols_p = 3;
    else
        n_rows_p = 3; n_cols_p = 3;
    end

    % Create figure for this priority
    priority_color_name = {'Red', 'Orange', 'Yellow', 'Yellow-green', ...
                           'Green', 'Blue', 'Purple', 'Violet'};

    fig_priority = figure('Name', sprintf('Priority %d Parameters', priority_level), ...
                          'Position', [150 + priority_level*20, 150 + priority_level*20, 1200, 800]);

    % Plot each parameter in this priority group
    for subplot_idx = 1:n_params_priority
        param_idx = priority_indices(subplot_idx);
        ptest = param_tests{param_idx};

        subplot(n_rows_p, n_cols_p, subplot_idx);
        hold on;

        n_values = length(ptest.values);
        colors_plot = lines(7); % Support up to 7 test values

        for j = 1:n_values
            plot(time, param_results{param_idx}(j, :), 'LineWidth', 2, ...
                 'Color', colors_plot(mod(j-1, 7)+1, :), ...
                 'DisplayName', sprintf('%.2f', ptest.values(j)));
        end

        xlabel('Time', 'FontSize', 10);
        ylabel('Fatigue (Fg)', 'FontSize', 10);
        title(ptest.label, 'FontSize', 11, 'Interpreter', 'tex');
        legend('Location', 'best', 'FontSize', 9);
        grid on;
        ylim([0, 1]);
        xlim([0, T]);
    end

    % Add overall title
    title_text = sprintf('Priority %d Parameters', priority_level);
    annotation(fig_priority, 'textbox', [0 0.96 1 0.04], ...
        'String', title_text, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% ===================== SAVE RESULTS =====================
save('sensitivity_parameters_complete.mat', 'param_results', 'param_tests', ...
     'sensitivity_summary', 'time', 'T', 'dt', 'params_baseline');

fprintf('✓ Results saved: sensitivity_parameters_complete.mat\n');
fprintf('✓ Visualizations created\n\n');

fprintf('=== ANALYSIS COMPLETE ===\n');
fprintf('Most sensitive parameter: %s (Range=%.4f)\n', sensitivity_summary(1).label, sensitivity_summary(1).range);
fprintf('Least sensitive parameter: %s (Range=%.4f)\n', sensitivity_summary(end).label, sensitivity_summary(end).range);
fprintf('\n');

end % end main function


%% ===================== HELPER FUNCTION =====================
function Fg = run_simulation(Te, Ti, Tc, Tp, Fi, Fl, Hr, Sh, St, numSteps, dt, params)
    % Extract parameters
    w_Td1 = params.w_Td1; w_Td2 = params.w_Td2; beta_Td = params.beta_Td;
    delta_Sq = params.delta_Sq;
    beta_Hl = params.beta_Hl; w_Hl1 = params.w_Hl1; w_Hl2 = params.w_Hl2;
    beta_Pp = params.beta_Pp; w_Pp1 = params.w_Pp1; w_Pp2 = params.w_Pp2; w_Pp3 = params.w_Pp3;
    alpha_Em = params.alpha_Em; eta_Em = params.eta_Em;
    w_Ge1 = params.w_Ge1; w_Ge2 = params.w_Ge2;
    w_Re1 = params.w_Re1; w_Re2 = params.w_Re2; w_Re3 = params.w_Re3;
    w_Sx1 = params.w_Sx1; w_Sx2 = params.w_Sx2; gamma_Sx = params.gamma_Sx;
    w_Sp1 = params.w_Sp1; w_Sp2 = params.w_Sp2; beta_Sp = params.beta_Sp;
    w_Rt1 = params.w_Rt1; w_Rt2 = params.w_Rt2;
    beta_Rt1 = params.beta_Rt1; beta_Rt2 = params.beta_Rt2; beta_Rt3 = params.beta_Rt3;

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
