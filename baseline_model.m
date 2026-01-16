clc;
clear;

%% ===================== TIME & SIM SETUP =====================
T  = 50;           % total duration (arbitrary units)
dt = 0.1;          % step size (0<dt<=1 is stable for these simple trackers)
numSteps  = T/dt;  % number of steps
time = linspace(0, T, numSteps);

%% ===================== INPUTS =====================
% Initialize arrays to store the variables over time (all normalized [0,1])
Te = zeros(1, numSteps); % Thermal environment (heat) a/
Ti = zeros(1, numSteps); % Training intensity a
Tc = zeros(1, numSteps); % Training complexity a
Tp = zeros(1, numSteps); % Training pressure b
Fi = zeros(1, numSteps); % Fluid intake b
Fl = zeros(1, numSteps); % Fitness level b
Hr = zeros(1, numSteps); % Heart-rate strain (0=low heart rate, 1=high heart rate) b a/
Sh = zeros(1, numSteps); % Sleep hours (normalized)b
St = zeros(1, numSteps); % Sleep timing / alignment (normalized)b

% Initial values
for t = 1:numSteps
    Te(t) = 0.9;     % Thermal environment
    Ti(t) = 0.9;     % Training intensity
    Tc(t) = 0.9;     % Training complexity
    Tp(t) = 0.9;     % Training pressure
    Hr(t) = 0.9;     % Heart rate
    Fi(t) = 0.8;     % Fluid intake
    Fl(t) = 0.2;     % Fitness level
    Sh(t) = 0.2;     % Sleep hours
    St(t) = 0.2;     % Sleep timing
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

% Recovery Effort (CALIBRATED for better moderate fatigue detection)
w_Re1 = 0.25;  % Reduced from 0.4 - less spare capacity recovery
w_Re2 = 0.20;  % Reduced from 0.3 - less effort-driven recovery
w_Re3 = 0.20;  % Reduced from 0.3 - less preparedness recovery

% Short-term Exhaustion (CALIBRATED to increase fatigue sensitivity)
w_Sx1 = 0.70;  % Reduced from 0.95（0.65） - less dependency on previous fatigue
w_Sx2 = 1.10;  % Increased from 0.95（1.15） - more sensitive to current effort
gamma_Sx = 0.50;  % Reduced from 0.8 - weaker recovery dampening

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

% NOTE: Initial values will be set after calculating Pp(1) for personalized baseline

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
% Adjust initial fatigue based on physical preparedness
initial_stress = 1 - Pp(1);  % Low preparedness → higher baseline stress
Fg(1) = 0.35 + 0.35 * initial_stress;  % Range: 0.35 (good prep) to 0.70 (poor prep)
Lx(1) = Fg(1);  % Match initial long-term exhaustion to fatigue
Lp(1) = 0.40 + 0.25 * initial_stress;  % Range: 0.40 to 0.65
Sx(1) = 0.40 + 0.30 * initial_stress;  % Range: 0.40 to 0.70

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

%% ===================== PLOTS =====================
minY = 0;
maxY = 1.05;

% 1) Short-term & Long-term stress/exhaustion + Fatigue
figure('Name','Temporal Relationships');
plot(time, Sp, 'LineWidth', 2);                    % Short-term Experienced Pressure
hold on;
plot(time, Lp, 'LineWidth', 2);                    % Long-term Pressure
plot(time, Sx, 'LineWidth', 2);                    % Short-term Exhaustion
plot(time, Lx, 'LineWidth', 2);                    % Long-term Exhaustion
plot(time, Fg, 'LineWidth', 2);                    % Fatigue
legend('Sp (Short Pressure)', ...
       'Lp (Long Pressure)', 'Sx (Short Exhaustion)','Lx (Long Exhaustion)','Fg (Fatigue)', ...
       'Location','best');
xlabel('Time'); ylabel('Level'); title('Short/Long-term Stress & Fatigue');
ylim([minY maxY]);
grid on;

% 2) Effort, Capacity, and Recovery
figure('Name','Effort & Capacity Dynamics');
plot(time, Em, 'LineWidth', 2);                    % Effort Motivation
hold on;
plot(time, Cp, 'LineWidth', 2);                    % Critical Point (capacity proxy)
plot(time, Ge, 'LineWidth', 2);                    % Generated Effort
plot(time, Re, 'LineWidth', 2);                    % Recovery Effort
plot(time, Rt, 'LineWidth', 2);                    % Recovert Time
legend('Em (Effort Motivation)','Cp (Critical Point)','Ge (Effort Output)','Re (Recovery Effort)', 'Rt (Recovery Time)', 'Location','best');
xlabel('Time'); ylabel('Level'); title('Effort, Capacity & Recovery');
ylim([minY maxY]);
grid on;

% 3) Inputs & Preparedness drivers
figure('Name','Task & Preparedness Drivers');
plot(time, Td, 'LineWidth', 2);                    % Task Demands
hold on;
plot(time, Hl, 'LineWidth', 2);                    % Hydration Level
plot(time, Pp, 'LineWidth', 2);                    % Physical Preparedness
plot(time, Sq, 'LineWidth', 2);                    % Sleep Quality
legend('Td (Task Demands)','Hl (Hydration Level)','Pp (Physical Preparedness)','Sq (Sleep Quality)', 'Location','best');
xlabel('Time'); ylabel('Level'); title('Task Demand & Physiological Drivers');
ylim([minY maxY]);
grid on;

% 4) Fatigue only (Fg)
figure('Name','Fatigue (Fg) Only');
plot(time, Fg, 'LineWidth', 2);
%legend('Fg (Fatigue)', 'Location','best');
xlabel('Time'); ylabel('Fatigue Level'); title('Fatigue Over Time (Fg)');
ylim([minY maxY]);
grid on;

