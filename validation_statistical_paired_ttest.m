% Statistical Validation - Paired t-Test Analysis
% ================================================
% This script performs paired t-test to compare simulated fatigue (Fg) with
% observed fatigue from questionnaire (FAS).
%
% Input:  validation_data.mat (from validate_with_questionnaire_data.m)
% Output: Console report, validation figure, and summary text file
%
% Author: Fatigue Simulation Research Team
% Date: 2025-12-26

function statistical_validation_paired_ttest()

clc;
fprintf('================================================================\n');
fprintf('STATISTICAL VALIDATION - PAIRED t-TEST ANALYSIS\n');
fprintf('================================================================\n\n');

%% ===================== LOAD VALIDATION DATA =====================
fprintf('Loading validation data...\n');

if ~exist('validation_data.mat', 'file')
    error('validation_data.mat not found! Please run validate_with_questionnaire_data.m first.');
end

load('validation_data.mat');

fprintf('[OK] Loaded data for %d respondents\n\n', n_respondents);

%% ===================== HYPOTHESIS DEVELOPMENT =====================
fprintf('================================================================\n');
fprintf('HYPOTHESIS DEVELOPMENT\n');
fprintf('================================================================\n\n');

fprintf('Null Hypothesis (H0):\n');
fprintf('  The mean difference between paired predicted and observed\n');
fprintf('  fatigue levels is zero.\n');
fprintf('  H0: mu_d = 0\n\n');

fprintf('Alternative Hypothesis (H1):\n');
fprintf('  The mean difference between paired predicted and observed\n');
fprintf('  fatigue levels is NOT zero.\n');
fprintf('  H1: mu_d ≠ 0\n\n');

fprintf('Test Type: Two-tailed paired t-test\n');
fprintf('Significance Level: alpha = 0.05\n\n');

%% ===================== DESCRIPTIVE STATISTICS =====================
fprintf('================================================================\n');
fprintf('DESCRIPTIVE STATISTICS\n');
fprintf('================================================================\n\n');

% Calculate differences
differences = Fg_simulated - Fatigue_observed;

% Calculate statistics
mean_diff = mean(differences);
std_diff = std(differences);
n = n_respondents;

fprintf('Sample Size (n): %d respondents\n\n', n);

fprintf('Simulated Fatigue (Fg_simulated):\n');
fprintf('  Mean:   %.4f\n', mean(Fg_simulated));
fprintf('  Std:    %.4f\n', std(Fg_simulated));
fprintf('  Range:  [%.4f, %.4f]\n\n', min(Fg_simulated), max(Fg_simulated));

fprintf('Observed Fatigue (Fatigue_observed from FAS):\n');
fprintf('  Mean:   %.4f\n', mean(Fatigue_observed));
fprintf('  Std:    %.4f\n', std(Fatigue_observed));
fprintf('  Range:  [%.4f, %.4f]\n\n', min(Fatigue_observed), max(Fatigue_observed));

fprintf('Differences (Simulated - Observed):\n');
fprintf('  Mean (mu_d):        %.4f\n', mean_diff);
fprintf('  Std Deviation (sd): %.4f\n', std_diff);
fprintf('  Range:              [%.4f, %.4f]\n\n', min(differences), max(differences));

%% ===================== PAIRED t-TEST =====================
fprintf('================================================================\n');
fprintf('PAIRED t-TEST CALCULATION\n');
fprintf('================================================================\n\n');

% Calculate t-statistic manually to show the formula
standard_error = std_diff / sqrt(n);
t_statistic = mean_diff / standard_error;
degrees_of_freedom = n - 1;

fprintf('Test Statistic Calculation:\n');
fprintf('  t = mu_d / (sd / sqrt(n))\n');
fprintf('  t = %.4f / (%.4f / sqrt(%d))\n', mean_diff, std_diff, n);
fprintf('  t = %.4f / %.4f\n', mean_diff, standard_error);
fprintf('  t = %.4f\n\n', t_statistic);

fprintf('Degrees of Freedom: df = n - 1 = %d - 1 = %d\n\n', n, degrees_of_freedom);

%% ===================== p-VALUE CALCULATION =====================
fprintf('================================================================\n');
fprintf('p-VALUE TESTING\n');
fprintf('================================================================\n\n');

% Calculate p-value manually using t-distribution
% For two-tailed test: p = 2 * P(T > |t|) where T ~ t(df)
% Use the cumulative distribution function (CDF)

% MATLAB's tcdf is in Statistics Toolbox, so we'll use a manual approximation
% or check if it's available
if exist('tcdf', 'file') == 2
    % Statistics Toolbox available
    p_value = 2 * (1 - tcdf(abs(t_statistic), degrees_of_freedom));
    h = (p_value <= 0.05);

    % Calculate confidence interval
    t_critical = tinv(0.975, degrees_of_freedom);  % 97.5th percentile for 95% CI
    margin = t_critical * standard_error;
    ci = [mean_diff - margin; mean_diff + margin];
else
    % Statistics Toolbox NOT available - use manual calculation
    % Approximate using normal distribution for df >= 10 (reasonable approximation)
    fprintf('Note: Statistics Toolbox not available. Using normal approximation.\n');
    fprintf('      (Valid for df >= 10, current df = %d)\n\n', degrees_of_freedom);

    % Standard normal approximation
    % P(Z > |z|) where Z ~ N(0,1)
    z_score = abs(t_statistic);

    % Approximate CDF using error function (erf is available in base MATLAB)
    % P(Z < z) = 0.5 * (1 + erf(z/sqrt(2)))
    p_one_tail = 0.5 * (1 - erf(z_score / sqrt(2)));
    p_value = 2 * p_one_tail;

    h = (p_value <= 0.05);

    % Approximate t-critical value (use z-critical for large df)
    z_critical = 1.96;  % For 95% CI (two-tailed)
    margin = z_critical * standard_error;
    ci = [mean_diff - margin; mean_diff + margin];
end

fprintf('Two-tailed p-value: %.6f\n\n', p_value);

fprintf('Interpretation:\n');
fprintf('  Significance level (alpha): 0.05\n');
fprintf('  p-value:                    %.6f\n', p_value);

if p_value > 0.05
    fprintf('  Result: p-value > 0.05\n\n');
else
    fprintf('  Result: p-value <= 0.05\n\n');
end

%% ===================== DECISION =====================
fprintf('================================================================\n');
fprintf('FINAL DECISION\n');
fprintf('================================================================\n\n');

if h == 0  % h=0 means fail to reject null hypothesis
    fprintf('DECISION: Fail to reject the null hypothesis (H0)\n\n');

    fprintf('CONCLUSION:\n');
    fprintf('  There is NO statistically significant difference between\n');
    fprintf('  the simulated fatigue levels (model predictions) and the\n');
    fprintf('  observed fatigue levels (questionnaire data).\n\n');

    fprintf('INTERPRETATION:\n');
    fprintf('  The fatigue simulation model produces predictions that are\n');
    fprintf('  statistically consistent with empirical questionnaire data\n');
    fprintf('  from %d Malaysian military cadets.\n\n', n);

    fprintf('  This validates that the model accurately captures real-world\n');
    fprintf('  fatigue patterns and can be trusted for predicting fatigue\n');
    fprintf('  in military training scenarios.\n\n');

    validation_status = 'VALIDATED';
else
    fprintf('DECISION: Reject the null hypothesis (H0)\n\n');

    fprintf('CONCLUSION:\n');
    fprintf('  There IS a statistically significant difference between\n');
    fprintf('  the simulated fatigue levels and observed fatigue levels.\n\n');

    fprintf('INTERPRETATION:\n');
    fprintf('  The model shows systematic bias (over-prediction or under-prediction).\n');
    fprintf('  Model calibration or parameter tuning may be needed to improve\n');
    fprintf('  agreement with empirical data.\n\n');

    if mean_diff > 0
        fprintf('  Bias Direction: Model OVER-predicts fatigue (Fg_sim > Fg_obs)\n\n');
    else
        fprintf('  Bias Direction: Model UNDER-predicts fatigue (Fg_sim < Fg_obs)\n\n');
    end

    validation_status = 'REQUIRES CALIBRATION';
end

fprintf('================================================================\n');
fprintf('VALIDATION STATUS: %s\n', validation_status);
fprintf('================================================================\n\n');

%% ===================== CONFIDENCE INTERVAL =====================
fprintf('95%% Confidence Interval for Mean Difference:\n');
fprintf('  [%.4f, %.4f]\n\n', ci(1), ci(2));

if ci(1) <= 0 && ci(2) >= 0
    fprintf('  The 95%% CI includes zero, supporting the conclusion that\n');
    fprintf('  there is no significant difference between simulated and observed fatigue.\n\n');
end

%% ===================== VISUALIZATION =====================
fprintf('Creating validation figure...\n');

figure('Name', 'Paired t-Test Validation', 'Position', [100, 100, 1200, 500]);

% Subplot 1: Observed vs Simulated Scatter Plot
subplot(1, 2, 1);
hold on;

% Plot perfect agreement line (y=x)
plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Perfect Agreement (y=x)');

% Plot data points
scatter(Fatigue_observed, Fg_simulated, 80, 'filled', 'MarkerFaceAlpha', 0.6, ...
        'DisplayName', 'Respondents (n=30)');

xlabel('Observed Fatigue (FAS Questionnaire)', 'FontSize', 11);
ylabel('Simulated Fatigue (Model Prediction)', 'FontSize', 11);
title('Observed vs Simulated Fatigue', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on;
axis([0 1 0 1]);
axis square;

% Add statistics text box
text_str = sprintf('n = %d\nMean diff = %.4f\nSD = %.4f\nt = %.4f\np = %.4f', ...
                   n, mean_diff, std_diff, t_statistic, p_value);
text(0.05, 0.92, text_str, 'FontSize', 9, 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black');

% Subplot 2: Differences Distribution
subplot(1, 2, 2);
hold on;

% Histogram of differences
histogram(differences, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'white', ...
          'FaceAlpha', 0.7, 'Normalization', 'probability');

% Add vertical line at zero
yl = ylim;
plot([0, 0], yl, 'r--', 'LineWidth', 2, 'DisplayName', 'Zero Difference');

% Add vertical line at mean difference
plot([mean_diff, mean_diff], yl, 'b-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Mean = %.4f', mean_diff));

xlabel('Difference (Simulated - Observed)', 'FontSize', 11);
ylabel('Probability', 'FontSize', 11);
title('Distribution of Differences', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;

% Add overall title
if exist('sgtitle', 'file') == 2
    sgtitle(sprintf('Statistical Validation: Paired t-Test | Status: %s', validation_status), ...
            'FontSize', 14, 'FontWeight', 'bold');
else
    annotation(gcf, 'textbox', [0 0.96 1 0.04], ...
        'String', sprintf('Statistical Validation: Paired t-Test | Status: %s', validation_status), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('[OK] Figure created!\n\n');

%% ===================== SAVE SUMMARY REPORT =====================
fprintf('Saving summary report...\n');

fid = fopen('validation_results.txt', 'w');

fprintf(fid, '================================================================\n');
fprintf(fid, 'STATISTICAL VALIDATION REPORT - PAIRED t-TEST\n');
fprintf(fid, '================================================================\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));

fprintf(fid, 'HYPOTHESIS:\n');
fprintf(fid, '  H0: mu_d = 0 (no difference between simulated and observed)\n');
fprintf(fid, '  H1: mu_d ≠ 0 (significant difference exists)\n\n');

fprintf(fid, 'SAMPLE:\n');
fprintf(fid, '  Sample size (n): %d respondents\n\n', n);

fprintf(fid, 'DESCRIPTIVE STATISTICS:\n');
fprintf(fid, '  Simulated Fatigue - Mean: %.4f, SD: %.4f\n', mean(Fg_simulated), std(Fg_simulated));
fprintf(fid, '  Observed Fatigue  - Mean: %.4f, SD: %.4f\n', mean(Fatigue_observed), std(Fatigue_observed));
fprintf(fid, '  Differences       - Mean: %.4f, SD: %.4f\n\n', mean_diff, std_diff);

fprintf(fid, 'TEST RESULTS:\n');
fprintf(fid, '  t-statistic:           %.4f\n', t_statistic);
fprintf(fid, '  Degrees of freedom:    %d\n', degrees_of_freedom);
fprintf(fid, '  p-value (two-tailed):  %.6f\n', p_value);
fprintf(fid, '  Significance level:    0.05\n\n');

fprintf(fid, '95%% CONFIDENCE INTERVAL:\n');
fprintf(fid, '  [%.4f, %.4f]\n\n', ci(1), ci(2));

fprintf(fid, 'DECISION:\n');
if h == 0
    fprintf(fid, '  Fail to reject H0\n\n');
else
    fprintf(fid, '  Reject H0\n\n');
end

fprintf(fid, 'VALIDATION STATUS: %s\n\n', validation_status);

fprintf(fid, 'INTERPRETATION:\n');
if h == 0
    fprintf(fid, '  The fatigue simulation model is statistically validated.\n');
    fprintf(fid, '  Model predictions are consistent with empirical questionnaire data.\n');
else
    fprintf(fid, '  The model shows systematic bias and requires calibration.\n');
    if mean_diff > 0
        fprintf(fid, '  Model tends to OVER-predict fatigue.\n');
    else
        fprintf(fid, '  Model tends to UNDER-predict fatigue.\n');
    end
end

fprintf(fid, '\n================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '================================================================\n');

fclose(fid);

fprintf('[OK] Summary report saved to: validation_results.txt\n\n');

%% ===================== SUMMARY =====================
fprintf('================================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('================================================================\n\n');

fprintf('Files generated:\n');
fprintf('  1. validation_results.txt - Statistical summary report\n');
fprintf('  2. Figure: Observed vs Simulated comparison\n\n');

fprintf('Key Results:\n');
fprintf('  Sample size:   %d\n', n);
fprintf('  Mean diff:     %.4f\n', mean_diff);
fprintf('  t-statistic:   %.4f\n', t_statistic);
fprintf('  p-value:       %.6f\n', p_value);
fprintf('  Status:        %s\n\n', validation_status);

fprintf('================================================================\n\n');

end % end main function
