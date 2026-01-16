% Questionnaire Data Validation - Simulation Runner
% ==================================================
% This script reads processed questionnaire data and runs the fatigue simulation
% model for each respondent to collect predicted fatigue values.
%
% Input:  Questionnaire_Data_Processed.csv (14 respondents with 9 normalized inputs)
% Output: validation_data.mat (simulated vs observed fatigue)
%
% Author: Fatigue Simulation Research Team
% Date: 2025-12-26

function validate_with_questionnaire_data()

clc;
fprintf('================================================================\n');
fprintf('FATIGUE MODEL VALIDATION WITH QUESTIONNAIRE DATA\n');
fprintf('================================================================\n\n');

%% ===================== READ PROCESSED QUESTIONNAIRE DATA =====================
csv_file = 'Questionnaire_Data_Processed.csv';

fprintf('Reading questionnaire data from: %s\n', csv_file);

% Read CSV file
data = readtable(csv_file);
n_respondents = height(data);

fprintf('[OK] Loaded %d respondents\n\n', n_respondents);

%% ===================== EXTRACT INPUT COLUMNS =====================
% Extract 9 normalized input variables (all in range [0,1])
Hr = data.Hr;  % Heart Rate
Sh = data.Sh;  % Sleep Hours
St = data.St;  % Sleep Timing
Tp = data.Tp;  % Training Pressure
Fi = data.Fi;  % Fluid Intake
Fl = data.Fl;  % Fitness Level
Te = data.Te;  % Thermal Environment
Ti = data.Ti;  % Training Intensity
Tc = data.Tc;  % Training Complexity

% Extract observed Fatigue from FAS questionnaire
Fatigue_observed = data.Fatigue;

% Extract metadata for reference
respondent_names = data.Name;
service_category = data.Service_Category;
level_of_training = data.Level_of_Training;

fprintf('Extracted inputs and observed fatigue for %d respondents\n\n', n_respondents);

%% ===================== SIMULATION PARAMETERS =====================
% Time setup (same as original_model.m)
T = 50;           % Total duration
dt = 0.1;         % Time step
numSteps = T/dt;  % Number of steps

% Model parameters (matching refined_model.m CALIBRATED version)
params = struct();

% Task Demands
params.w_Td1 = 0.4;
params.w_Td2 = 0.6;
params.beta_Td = 0.7;

% Sleep Quality
params.delta_Sq = 0.7;

% Hydration Level
params.beta_Hl = 0.5;
params.w_Hl1 = 0.6;
params.w_Hl2 = 0.4;

% Physical Preparedness
params.beta_Pp = 0.7;
params.w_Pp1 = 0.25;
params.w_Pp2 = 0.30;
params.w_Pp3 = 0.45;

% Effort Motivation
params.alpha_Em = 0.4;
params.eta_Em = 0.5;

% Generated Effort
params.w_Ge1 = 0.85;
params.w_Ge2 = 0.15;

% Recovery Effort (CALIBRATED)
params.w_Re1 = 0.25;
params.w_Re2 = 0.20;
params.w_Re3 = 0.20;

% Short-term Exhaustion (CALIBRATED)
params.w_Sx1 = 0.70;
params.w_Sx2 = 1.10;
params.gamma_Sx = 0.50;

% Short-term Experienced Pressure
params.w_Sp1 = 0.6;
params.w_Sp2 = 0.4;
params.beta_Sp = 0.7;

% Recovery Time
params.w_Rt1 = 0.5;
params.w_Rt2 = 0.5;
params.beta_Rt1 = 0.3;
params.beta_Rt2 = 0.3;
params.beta_Rt3 = 0.4;

%% ===================== RUN SIMULATIONS FOR ALL RESPONDENTS =====================
fprintf('Running fatigue simulation for each respondent...\n');
fprintf('Simulation settings: T=%.1f, dt=%.2f, numSteps=%d\n\n', T, dt, numSteps);

Fg_simulated = zeros(n_respondents, 1);

for i = 1:n_respondents
    fprintf('  [%2d/%2d] Simulating: %s... ', i, n_respondents, respondent_names{i});

    % Create constant input arrays for this respondent
    Te_input = ones(1, numSteps) * Te(i);
    Ti_input = ones(1, numSteps) * Ti(i);
    Tc_input = ones(1, numSteps) * Tc(i);
    Tp_input = ones(1, numSteps) * Tp(i);
    Fi_input = ones(1, numSteps) * Fi(i);
    Fl_input = ones(1, numSteps) * Fl(i);
    Hr_input = ones(1, numSteps) * Hr(i);
    Sh_input = ones(1, numSteps) * Sh(i);
    St_input = ones(1, numSteps) * St(i);

    % Run simulation
    Fg_output = run_fatigue_simulation(Te_input, Ti_input, Tc_input, Tp_input, ...
                                       Fi_input, Fl_input, Hr_input, Sh_input, St_input, ...
                                       numSteps, dt, params);

    % Store final fatigue value
    Fg_simulated(i) = Fg_output(end);

    fprintf('Fg_sim=%.4f, Fg_obs=%.4f\n', Fg_simulated(i), Fatigue_observed(i));
end

fprintf('\n[OK] All simulations completed!\n\n');

%% ===================== DETAILED COMPARISON TABLE =====================
fprintf('================================================================\n');
fprintf('DETAILED COMPARISON: OBSERVED vs SIMULATED FATIGUE\n');
fprintf('================================================================\n\n');

fprintf('%-40s %12s %12s %12s\n', 'Respondent Name', 'Observed', 'Simulated', 'Difference');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:n_respondents
    diff = Fg_simulated(i) - Fatigue_observed(i);
    fprintf('%-40s %12.4f %12.4f %12.4f\n', ...
            respondent_names{i}, Fatigue_observed(i), Fg_simulated(i), diff);
end

fprintf('--------------------------------------------------------------------------------\n');

% Calculate statistics for the detailed table
mean_observed = mean(Fatigue_observed);
mean_simulated = mean(Fg_simulated);
mean_diff = mean(Fg_simulated - Fatigue_observed);

std_observed = std(Fatigue_observed);
std_simulated = std(Fg_simulated);
std_diff = std(Fg_simulated - Fatigue_observed);

fprintf('%-40s %12.4f %12.4f %12.4f\n', 'MEAN', mean_observed, mean_simulated, mean_diff);
fprintf('%-40s %12.4f %12.4f %12.4f\n', 'STD', std_observed, std_simulated, std_diff);

fprintf('\n');

%% ===================== CALCULATE DIFFERENCES =====================
differences = Fg_simulated - Fatigue_observed;

fprintf('================================================================\n');
fprintf('VALIDATION DATA SUMMARY\n');
fprintf('================================================================\n\n');

fprintf('Respondents: %d\n\n', n_respondents);

fprintf('%-30s %10s %10s %10s %10s\n', 'Variable', 'Min', 'Max', 'Mean', 'Std');
fprintf('------------------------------------------------------------------------\n');
fprintf('%-30s %10.4f %10.4f %10.4f %10.4f\n', 'Fg_simulated', min(Fg_simulated), max(Fg_simulated), mean(Fg_simulated), std(Fg_simulated));
fprintf('%-30s %10.4f %10.4f %10.4f %10.4f\n', 'Fatigue_observed', min(Fatigue_observed), max(Fatigue_observed), mean(Fatigue_observed), std(Fatigue_observed));
fprintf('%-30s %10.4f %10.4f %10.4f %10.4f\n', 'Differences (Sim - Obs)', min(differences), max(differences), mean(differences), std(differences));

fprintf('\n================================================================\n\n');

%% ===================== SAVE RESULTS =====================
save('validation_data.mat', 'Fg_simulated', 'Fatigue_observed', 'differences', ...
     'respondent_names', 'service_category', 'level_of_training', ...
     'n_respondents', 'T', 'dt', 'params');

fprintf('[OK] Results saved to: validation_data.mat\n');
fprintf('[OK] Ready for statistical testing!\n\n');

fprintf('Next step: Run validation_statistical_paired_ttest.m\n\n');

end % end main function


%% ===================== HELPER FUNCTION: FATIGUE SIMULATION =====================
function Fg = run_fatigue_simulation(Te, Ti, Tc, Tp, Fi, Fl, Hr, Sh, St, numSteps, dt, params)
    % Run fatigue simulation with given inputs and parameters.
    % Returns the Fg (Fatigue) time series.

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

    % Step 1 (t=1) instantaneous compute
    Td(1) = beta_Td * Ti(1) * (w_Td1*Tc(1) + w_Td2*Tp(1)) + (1 - beta_Td) * (Te(1) * Ti(1));
    Sq(1) = delta_Sq * Sh(1) + (1-delta_Sq) * St(1);
    Hl(1) = (beta_Hl * Fi(1)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(1) + w_Hl2 * Td(1)));
    Pp(1) = beta_Pp * (w_Pp1*Hl(1) + w_Pp2*Fl(1) + w_Pp3*Sq(1)) + (1-beta_Pp) * (1 - Hr(1));

    % CALIBRATED INITIAL CONDITIONS (personalized based on preparedness)
    initial_stress = 1 - Pp(1);
    Fg(1) = 0.35 + 0.35 * initial_stress;
    Lx(1) = Fg(1);
    Lp(1) = 0.40 + 0.25 * initial_stress;
    Sx(1) = 0.40 + 0.30 * initial_stress;

    Cp(1) = Pp(1) * (1 - Sx(1));
    Em(1) = exp( -eta_Em * ((Td(1) - (1-Lp(1)) ) / alpha_Em) ^ 2 );
    Ge(1) = Em(1) * (w_Ge1*Td(1) + w_Ge2*Cp(1));
    Re(1) = w_Re1 * min(1, max(0, Cp(1) - Ge(1))) + w_Re2 * Ge(1) + w_Re3 * min(1, max(0, (Pp(1) - Cp(1)) / max(Pp(1), 1e-9)));
    Sx(1) = (w_Sx1*Fg(1) + w_Sx2*max(0, Ge(1)-Cp(1))) * (1 - gamma_Sx*Re(1));
    Sp(1) = beta_Sp * (w_Sp1 * Sx(1) + w_Sp2 * max(0, Ge(1)-Cp(1))) + (1-beta_Sp) * (1-Pp(1));
    Rt(1) = (w_Rt1*Sx(1) + w_Rt2*Fg(1)) * (1 - (beta_Rt1*Pp(1) + beta_Rt2*Sq(1) + beta_Rt3*Re(1)));

    % Main loop (t=2 to numSteps)
    for t = 2:numSteps
        % Instantaneous layer
        Td(t) = beta_Td * Ti(t) * (w_Td1*Tc(t) + w_Td2*Tp(t)) + (1 - beta_Td) * (Te(t) * Ti(t));
        Sq(t) = delta_Sq * Sh(t) + (1-delta_Sq) * St(t-1);
        Hl(t) = (beta_Hl * Fi(t)) + (1 - beta_Hl) * (1 - (w_Hl1 * Te(t) + w_Hl2 * Td(t)));
        Pp(t) = beta_Pp * (w_Pp1*Hl(t) + w_Pp2*Fl(t) + w_Pp3*Sq(t)) + (1-beta_Pp) * (1 - Hr(t));
        Cp(t) = Pp(t) * (1 - Sx(t-1));
        Em(t) = exp( -eta_Em * ((Td(t) - (1-Lp(t-1)) ) / alpha_Em) ^ 2 );
        Ge(t) = Em(t) * (w_Ge1*Td(t) + w_Ge2*Cp(t));
        Re(t) = w_Re1 * min(1, max(0, Cp(t) - Ge(t))) + w_Re2 * Ge(t) + w_Re3 * min(1, max(0, (Pp(t) - Cp(t)) / max(Pp(t), 1e-9)));
        Sx(t) = (w_Sx1*Fg(t-1) + w_Sx2*max(0, Ge(t)-Cp(t))) * (1 - gamma_Sx*Re(t));
        Sp(t) = beta_Sp * (w_Sp1 * Sx(t-1) + w_Sp2 * max(0, Ge(t)-Cp(t))) + (1-beta_Sp) * (1-Pp(t));

        % Temporal layer (logistic growth trackers)
        Lp(t) = Lp(t-1) + (Sp(t) - Lp(t-1)) * Lp(t-1) * (1-Lp(t-1)) * dt;
        Lx(t) = Lx(t-1) + (Sx(t-1) - Lx(t-1)) * Lx(t-1) * (1-Lx(t-1)) * dt;
        Fg(t) = Fg(t-1) + (Lx(t-1) - Fg(t-1)) * Fg(t-1) * (1-Fg(t-1)) * dt;

        Rt(t) = (w_Rt1*Sx(t-1) + w_Rt2*Fg(t-1)) * (1 - (beta_Rt1*Pp(t) + beta_Rt2*Sq(t) + beta_Rt3*Re(t)));
    end
end
