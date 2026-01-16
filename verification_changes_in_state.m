function fatigue_changes_in_state_cases()
    clc; close all; clear;

    % ===================== TIME SETUP =====================
    T  = 500;
    dt = 0.1;
    time = 0:dt:T;                 % 0,0.1,...,500  (5001 points)
    numSteps = numel(time);

    % ===================== PARAMETERS (same as your code) =====================
    P = getParams();

    % ===================== CASE A: Fg -> 1 =====================
    caseA.name = 'Case A (Fg -> 1): High load, low recovery';
    caseA.inputs.Te = 1; caseA.inputs.Ti = 1; caseA.inputs.Tc = 1; caseA.inputs.Tp = 1;
    caseA.inputs.Hr = 1; caseA.inputs.Fi = 0; caseA.inputs.Fl = 0; caseA.inputs.Sh = 0; caseA.inputs.St = 0;

    % Initial states (interior)
    caseA.init.Lp = 0.50;
    caseA.init.Lx = 0.50;
    caseA.init.Fg = 0.50;
    caseA.init.Sx = 0.50;

    outA = simulateFatigue(time, dt, numSteps, P, caseA.inputs, caseA.init);
    plotCase(time, outA, caseA.name);

    % ===================== CASE B: Fg -> 0 =====================
    caseB.name = 'Case B (Fg -> 0): Low load, strong recovery';
    % Make Task Demands ~ 0 and Preparedness ~ 1
    caseB.inputs.Te = 0; caseB.inputs.Ti = 0; caseB.inputs.Tc = 0; caseB.inputs.Tp = 0;
    caseB.inputs.Hr = 0; caseB.inputs.Fi = 1; caseB.inputs.Fl = 1; caseB.inputs.Sh = 1; caseB.inputs.St = 1;

    % Lock Lp at boundary 0 (so Em becomes small when Td far from 1)
    caseB.init.Lp = 1e-3;
    caseB.init.Lx = 0.10;
    caseB.init.Fg = 0.10;
    caseB.init.Sx = 0.10;

    outB = simulateFatigue(time, dt, numSteps, P, caseB.inputs, caseB.init);
    plotCase(time, outB, caseB.name);

    fprintf('\nDone. Two figures created (Fg->1 and Fg->0).\n');
end


% =====================================================================
% SIMULATION (same structure as your model, but with constant inputs)
% =====================================================================
function out = simulateFatigue(time, dt, numSteps, P, U, init)
    % Inputs as arrays (constant over time)
    Te = U.Te * ones(1,numSteps);
    Ti = U.Ti * ones(1,numSteps);
    Tc = U.Tc * ones(1,numSteps);
    Tp = U.Tp * ones(1,numSteps);
    Hr = U.Hr * ones(1,numSteps);
    Fi = U.Fi * ones(1,numSteps);
    Fl = U.Fl * ones(1,numSteps);
    Sh = U.Sh * ones(1,numSteps);
    St = U.St * ones(1,numSteps);

    % States
    Td = zeros(1,numSteps);
    Hl = zeros(1,numSteps);
    Sq = zeros(1,numSteps);
    Pp = zeros(1,numSteps);
    Em = zeros(1,numSteps);
    Cp = zeros(1,numSteps);
    Ge = zeros(1,numSteps);
    Re = zeros(1,numSteps);
    Sp = zeros(1,numSteps);
    Sx = zeros(1,numSteps);
    Lp = zeros(1,numSteps);
    Lx = zeros(1,numSteps);
    Fg = zeros(1,numSteps);

    % Init
    Lp(1) = init.Lp;
    Lx(1) = init.Lx;
    Fg(1) = init.Fg;
    Sx(1) = init.Sx;

    % Instant layer at t=1 (use same structure as your code)
    Td(1) = P.beta_Td * Ti(1) * (P.w_Td1*Tc(1) + P.w_Td2*Tp(1)) + (1 - P.beta_Td) * (Te(1) * Ti(1));
    Sq(1) = P.delta_Sq*Sh(1) + (1-P.delta_Sq)*St(1);
    Hl(1) = (P.beta_Hl*Fi(1)) + (1-P.beta_Hl)*(1 - (P.w_Hl1*Te(1) + P.w_Hl2*Td(1)));
    Pp(1) = P.beta_Pp*(P.w_Pp1*Hl(1) + P.w_Pp2*Fl(1) + P.w_Pp3*Sq(1)) + (1-P.beta_Pp)*(1-Hr(1));
    Cp(1) = Pp(1)*(1 - Sx(1));
    Em(1) = exp( -P.eta_Em * ((Td(1) - (1-Lp(1))) / P.alpha_Em)^2 );
    Ge(1) = Em(1)*(P.w_Ge1*Td(1) + P.w_Ge2*Cp(1));

    term1 = min(1, max(0, Cp(1)-Ge(1)));
    term3 = min(1, max(0, (Pp(1)-Cp(1)) / max(Pp(1),1e-9)));
    Re(1) = P.w_Re1*term1 + P.w_Re2*Ge(1) + P.w_Re3*term3;

    overload1 = max(0, Ge(1)-Cp(1));
    Sx(1) = (P.w_Sx1*Fg(1) + P.w_Sx2*overload1)*(1 - P.gamma_Sx*Re(1));
    Sp(1) = P.beta_Sp*(P.w_Sp1*Sx(1) + P.w_Sp2*overload1) + (1-P.beta_Sp)*(1-Pp(1));

    % Main loop
    for t = 2:numSteps
        Td(t) = P.beta_Td * Ti(t) * (P.w_Td1*Tc(t) + P.w_Td2*Tp(t)) + (1 - P.beta_Td) * (Te(t) * Ti(t));
        Sq(t) = P.delta_Sq*Sh(t) + (1-P.delta_Sq)*St(t);
        Hl(t) = (P.beta_Hl*Fi(t)) + (1-P.beta_Hl)*(1 - (P.w_Hl1*Te(t) + P.w_Hl2*Td(t)));
        Pp(t) = P.beta_Pp*(P.w_Pp1*Hl(t) + P.w_Pp2*Fl(t) + P.w_Pp3*Sq(t)) + (1-P.beta_Pp)*(1-Hr(t));

        Cp(t) = Pp(t)*(1 - Sx(t-1));

        % NOTE: Lp appears here (this is why long-term affects instantaneous)
        Em(t) = exp( -P.eta_Em * ((Td(t) - (1-Lp(t-1))) / P.alpha_Em)^2 );

        Ge(t) = Em(t)*(P.w_Ge1*Td(t) + P.w_Ge2*Cp(t));

        term1 = min(1, max(0, Cp(t)-Ge(t)));
        term3 = min(1, max(0, (Pp(t)-Cp(t)) / max(Pp(t),1e-9)));
        Re(t) = P.w_Re1*term1 + P.w_Re2*Ge(t) + P.w_Re3*term3;

        overload = max(0, Ge(t)-Cp(t));
        Sx(t) = (P.w_Sx1*Fg(t-1) + P.w_Sx2*overload)*(1 - P.gamma_Sx*Re(t));
        Sp(t) = P.beta_Sp*(P.w_Sp1*Sx(t-1) + P.w_Sp2*overload) + (1-P.beta_Sp)*(1-Pp(t));

        % Long-term trackers
        Lp(t) = Lp(t-1) + (Sp(t)-Lp(t-1))*Lp(t-1)*(1-Lp(t-1))*dt;
        Lx(t) = Lx(t-1) + (Sx(t-1)-Lx(t-1))*Lx(t-1)*(1-Lx(t-1))*dt;
        Fg(t) = Fg(t-1) + (Lx(t-1)-Fg(t-1))*Fg(t-1)*(1-Fg(t-1))*dt;
    end

    out.time = time;
    out.inputs = U;
    out.Td = Td; out.Pp = Pp; out.Em = Em; out.Ge = Ge; out.Re = Re;
    out.Sx = Sx; out.Sp = Sp; out.Lp = Lp; out.Lx = Lx; out.Fg = Fg;

    % "Changes in state" signals (like case study "changes in long-term stress")
    out.dFg = [0 diff(Fg)];
    out.dLx = [0 diff(Lx)];
end


% =====================================================================
% PLOTTING (case-study style: stimuli + change + state comparison)
% =====================================================================
function plotCase(time, out, figTitle)
    figure('Name', figTitle);

    % ---------- prep for sign check ----------
    gap = out.Lx - out.Fg;                  % Lx(t) - Fg(t)
    gap_prev = [gap(1) gap(1:end-1)];       % matches ΔFg(t) uses (t-1)

    % Sign consistency score (ignore exact zeros)
    s1 = sign(out.dFg);
    s2 = sign(gap_prev);
    mask = (s1 ~= 0) & (s2 ~= 0);
    if any(mask)
        agree_pct = 100 * mean(s1(mask) == s2(mask));
    else
        agree_pct = NaN;
    end

    % ===================== (1) Inputs =====================
    subplot(2,2,1);
    plot(time, out.inputs.Te*ones(size(time)), 'LineWidth', 2); hold on;
    plot(time, out.inputs.Ti*ones(size(time)), 'LineWidth', 2);
    plot(time, out.inputs.Sh*ones(size(time)), 'LineWidth', 2);
    plot(time, out.inputs.Fi*ones(size(time)), 'LineWidth', 2);
    title('Stimuli / Inputs (constant)');
    xlabel('Time'); ylabel('Level'); grid on; ylim([-0.05 1.05]);
    legend({'Te','Ti','Sh','Fi'}, 'Location','northeast');

    % ===================== (2) ΔFg =====================
    subplot(2,2,2);
    plot(time, out.dFg, 'LineWidth', 2);
    title('\DeltaFg(t) = Fg(t)-Fg(t-1)');
    xlabel('Time'); ylabel('\DeltaFg'); grid on;

    % ===================== (3) Lx vs Fg =====================
    subplot(2,2,3);
    plot(time, out.Lx, '--', 'LineWidth', 2); hold on;
    plot(time, out.Fg,  '-', 'LineWidth', 2);
    title('Driver vs State: Lx (dashed) and Fg (solid)');
    xlabel('Time'); ylabel('Level'); grid on; ylim([-0.05 1.05]);
    legend({'Lx (driver)','Fg (state)'}, 'Location','northeast');

    % ===================== (4) Chain =====================
    subplot(2,2,4);
    plot(time, out.Sx, 'LineWidth', 2); hold on;
    plot(time, out.Lx, 'LineWidth', 2);
    plot(time, out.Fg, 'LineWidth', 2);
    title('Chain: Sx, Lx, Fg');
    xlabel('Time'); ylabel('Level'); grid on; ylim([-0.05 1.05]);
    legend({'Sx','Lx','Fg'}, 'Location','northeast');

    sgtitle_compat(figTitle);
end


% =====================================================================
% sgtitle compatibility for Octave/MATLAB
% =====================================================================
function sgtitle_compat(str)
    if exist('sgtitle','file') == 2
        sgtitle(str);
        return;
    end

    if exist('suptitle','file') == 2
        suptitle(str);
        return;
    end

    % Pure Octave fallback (works without extra packages)
    annotation('textbox', [0 0.965 1 0.035], ...
        'String', str, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
end



% =====================================================================
% PARAMETERS (copied from your model)
% =====================================================================
function P = getParams()
    % Task Demands
    P.w_Td1 = 0.4; P.w_Td2 = 0.6; P.beta_Td = 0.7;

    % Sleep Quality
    P.delta_Sq = 0.7;

    % Hydration Level
    P.beta_Hl = 0.5; P.w_Hl1 = 0.6; P.w_Hl2 = 0.4;

    % Physical Preparedness
    P.beta_Pp = 0.7; P.w_Pp1 = 0.25; P.w_Pp2 = 0.30; P.w_Pp3 = 0.45;

    % Effort Motivation
    P.alpha_Em = 0.4; P.eta_Em = 0.5;

    % Generated Effort
    P.w_Ge1 = 0.85; P.w_Ge2 = 0.15;

    % Recovery Effort (CALIBRATED)
    P.w_Re1 = 0.25; P.w_Re2 = 0.20; P.w_Re3 = 0.20;

    % Short-term Exhaustion (CALIBRATED)
    P.w_Sx1 = 0.70; P.w_Sx2 = 1.10; P.gamma_Sx = 0.50;

    % Short-term Experienced Pressure
    P.w_Sp1 = 0.6; P.w_Sp2 = 0.4; P.beta_Sp = 0.7;
end


