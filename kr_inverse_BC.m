%% =========================================================================
%  HISTORY MATCHING USING CHIERICI RELATIVE PERMEABILITY MODEL
%  =========================================================================
%  Description:
%    This script estimates Chierici model parameters (A, L, B, M)
%    by matching fractional flow and pressure drop data generated from a
%    high‑fidelity analytical layered model.
%
%  Workflow:
%    1. Load historical data (PVI, fw_history, dp_history).
%    2. Define the Chierici parameter bounds and initial guess.
%    3. For each candidate parameter set:
%       a. Compute fractional flow and derivative using relativePermeability.
%       b. Find shock front saturation (Welge tangent) using findShockFront.
%       c. Compute breakthrough PVI = 1 / f'(sw_front).
%       d. For each historical PVI, compute outlet saturation (Buckley‑Leverett),
%          then fractional flow, and pressure drop using computePressureDrop.
%       e. Compute objective = sum( ( (fw_pred-fw_hist)/max(fw_hist) )^2 ) +
%                              sum( ( (dp_pred-dp_hist)/max(dp_hist) )^2 ).
%    4. Run optimization (fmincon) to minimize the objective.
%    5. Plot the final match.
%
% INPUTS
% This code takes inputs as 2 patches: 1 and2
% Patch 1: contains analytical results data
% Patch 2: initial guess and bounds of matching parameters
%
%  Author:   AbdAllah AbdRabElNaby Youssef ElSaieed
%  Date:     14th May 2026
%% =========================================================================

clear; clc; close all; warning off

% =========================================================================
%% PARSE INPUTS
% =========================================================================

% --- INPUTS: Patch 1 (Analytical results data) ---
%  Description:
%    The code load the analytical results saved after running
%    "MultiLayer.m". The results are saved in structure "combined" in
%    addtion to 'mu_nw', 'mu_w'.
%    analytical_results_file: .mat file that saves results of MultiLayer.m.
%    The user is reuired to change the name based on his/her naming 
% SYNTAX:
%    - analytical_results_file       :  .mat file that contains analytical results

analytical_results_file = 'Otway_miscible_results.mat';      % replace with .mat data produced from MultiLayer.m


% --- INPUTS: Patch 2 (Initial guess and bounds of matching parameters) ---
% SYNTAX:
%    -lb       : lower bound of matching parameters
%    -ub      : upper bound of matching parameters
%    -initial_guess      : initial guess of matching parameters
%    -injected_pv_limit      : set limit of matching region to exclude very large pore 

% Lower bounds for parameters [A, L, B, M]
lb = 0.0;
% Upper bounds
ub = 10.0;

% Initial guess (reasonable starting point)
initial_guess =  3.1547;

% Maximum value of injected pore volume to set match
injected_pv_limit = 2.2;



% =========================================================================
%% LOADING INPUTS
% =========================================================================
load(analytical_results_file,  'combined', 'mu_nw', 'mu_w', 'Sw_tanget_start', 'layers',...
    'Process', 'Miscible', 'Snwr', 'Swr', 'Krnw_end', 'Krw_end', 'SwInitial');

% Get layer parameters
porosity =  [layers.porosity]';
thickness = [layers.thickness]';
perm_abs = [layers.perm_abs]';


% Extract IPV, pressure drop and Fw
injected_pv_analytical = combined.injected_pv;        % PVI estimated analytically
pressure_drop_analytical_input = combined.pressure_drop;         % pressure drop esimated analytically 
F_w_analytical = combined.F_w;          % fractional flow estiamted analytically

% Ensure all are column vectors
injected_pv_analytical = injected_pv_analytical(:);
pressure_drop_analytical_input   = pressure_drop_analytical_input(:);
F_w_analytical   = F_w_analytical(:);
Sw_tanget_start = Sw_tanget_start(:);
Snwr = Snwr(:);
Swr = Swr(:);
Krnw_end = Krnw_end(:);
Krw_end = Krw_end(:);
SwInitial = SwInitial(:);

% =========================================================================
%% TRIMMING DATA AND EVALUATING ENDPOINTS
% =========================================================================

% Trimming data to certain region which is important
pressure_drop_analytical_input(injected_pv_analytical > injected_pv_limit) = [];
F_w_analytical(injected_pv_analytical > injected_pv_limit) = [];
injected_pv_analytical(injected_pv_analytical > injected_pv_limit) = [];

% Number of data points
n_points = length(injected_pv_analytical);

% Estimating effective permeability and endpoints
effective.perm_abs = sum(perm_abs .* thickness) / sum(thickness);    % effective absolute permeability
effective.thickness = sum(thickness);                                % net pay thickness
effective.Swr = sum(Swr .* porosity .* thickness) / sum(porosity .* thickness);   % effective wetting residual saturation
effective.Snwr = sum(Snwr .* porosity .* thickness) / sum(porosity .* thickness);   % effective non-wetting residual saturation
effective.Krnw_end = sum(Krnw_end .* perm_abs .* thickness) / sum(perm_abs .* thickness);   % effective non-wetting endpoint relative permeability
effective.Krw_end = sum(Krw_end .* perm_abs .* thickness) / sum(perm_abs .* thickness);   % effective wetting endpoint relative permeability
effective.SwInitial = sum(SwInitial .* porosity .* thickness) / sum(porosity .* thickness);   % effective initial saturation

% Determine start of tanget 
effective.Sw_tanget_start = sum(Sw_tanget_start .* thickness .* porosity) / (sum(thickness .* porosity));

% Determine conditions of flooding
effective.Miscible = Miscible;
effective.Process = Process;


% =========================================================================
%% CONSTRUCT OBJECTIVE FUNCTION BY DEFINING FORWARD MODEL (HOMOGENIZED SIMULATION)
% =========================================================================
% In this step pressure drop and fractional flow at the same PVI points as the
% analytical solution are evalauted for homogeneous model
% Anonymous function that passes extra parameters to objective
obj_fun = @(theta) objective(theta, injected_pv_analytical, pressure_drop_analytical_input, F_w_analytical, mu_w, mu_nw, effective);


% =========================================================================
%% RUN OPTIMIZATION
% =========================================================================
% In this step the optimium kr parameters are evaluated
% Use constrained nonlinear optimization (fmincon) with SQP algorithm.

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...                 % show iteration details
    'Algorithm', 'sqp', ...                % sequential quadratic programming
    'MaxFunctionEvaluations', 5000, ...    % maximum number of objective evaluations
    'StepTolerance', 1e-10, ...             % tolerance on parameter change
    'OptimalityTolerance', 1e-10);          % tolerance on first-order optimality



% Run fmincon
[theta_opt, J_opt, exitflag, output] = fmincon(obj_fun, initial_guess, [], [], [], [], lb, ub , [], options);

%% 7. DISPLAY OPTIMIZED PARAMETERS
% -------------------------------------------------------------------------
fprintf('\n========== OPTIMIZATION RESULTS ==========\n');
fprintf('Objective function value: %.6e\n', J_opt);
fprintf('Exit flag: %d  (%s)\n', exitflag, output.message);
fprintf('\nOptimized BC parameter:\n');
fprintf('  \lamda  = %.4f\n', theta_opt);

%% 8. EVALUATE FINAL MATCH AND PLOT RESULTS
% -------------------------------------------------------------------------
% Run forward model one last time with optimized parameters
[dp_opt, fw_opt] = forwardHomogenized(theta_opt, injected_pv_analytical, mu_w, mu_nw, effective);

% % Compute final relative errors
% eps_dp_final = (dp_opt - pressure_drop_analytical_input) ./ (pressure_drop_analytical_input + eps);
% eps_fw_final = (fw_opt - F_w_analytical) ./ (F_w_analytical + eps);
% max_rel_error_dp = max(abs(eps_dp_final));
% max_rel_error_fw = max(abs(eps_fw_final));
% 
% fprintf('\nFinal relative errors:\n');
% fprintf('  Max |ε_dp| = %.4f\n', max_rel_error_dp);
% fprintf('  Max |ε_fw| = %.4f\n', max_rel_error_fw);

% Create comparison plots
figure('Name', 'Inverse Modeling Results', 'Position', [100, 100, 1200, 500]);

subplot(1,2,1);
plot(injected_pv_analytical, F_w_analytical, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(injected_pv_analytical, fw_opt, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('PVI (Pore Volumes Injected)');
ylabel('Fractional Flow (f_w)');
title('Fractional Flow Match');
legend('Analytical (layered)', 'Homogenized (optimized)', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(injected_pv_analytical, pressure_drop_analytical_input, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(injected_pv_analytical, dp_opt, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('PVI (Pore Volumes Injected)');
ylabel('Pressure Drop (dimensionless or absolute)');
title('Pressure Drop Match');
legend('Analytical (layered)', 'Homogenized (optimized)', 'Location', 'best');
grid on;

sgtitle('Comparison between Analytical Layered Model and Optimized Homogenized Model');

%% 9. (OPTIONAL) SAVE RESULTS
% -------------------------------------------------------------------------

% save('inverse_modeling_results.mat', 'theta_opt', 'J_opt', 'dp_opt', 'fw_opt');


% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [dp_effective, fw_effective] = forwardHomogenized(theta, PVI_analytical, mu_w, mu_nw, effective)
% FORWARDHOMOGENIZED Pressure drop and outlet fractional flow of the
% equivalent homogeneous model at the PVI values of the analytical solution.
   
    
    % Unpack parameters and build relative permeability and fractional flow functions
    params.Swr = effective.Swr;
    params.Snwr = effective.Snwr;
    params.Krw_end = effective.Krw_end;
    params.Krnw_end = effective.Krnw_end;
    params.lambda = theta;
    
  

    % Run forward model to get fractional flow and total mobility
    [krw_func, krnw_func, ~, ~, mob_t_func, fw_func, dfw_dsw_func] = ...
      relativePermeability('BC', params, mu_w, mu_nw);
    
    % Find shock front saturation
    [sw_front, ~, slope] = findShockFront(effective.Swr, effective.Snwr, fw_func, dfw_dsw_func, ...
    effective.Sw_tanget_start, 'Process', effective.Process, 'Miscible', effective.Miscible, 'Verbose', false, 'Plot', false);
    
    % Determine saturation bounds based on process
    if strcmpi(effective.Process, 'drainage')
        sw_max = sw_front;  % For drainage, maximum is at front
        sw_min = effective.Swr;
    else
        sw_max = 1 - effective.Snwr;  % For imbibition, maximum is 1-Snwr
        sw_min = sw_front;
    end

    % Compute gamma parameter before breakthrough
    if effective.Miscible
        [gamma_breakthrough, ~, ~] = computeGammaParameter3(effective.Swr, effective.Snwr, sw_front, [], ...
            fw_func, dfw_dsw_func, mob_t_func, 'Verbose', false, 'Process', effective.Process, 'Compute', 'miscible', 'PlotIntegrand', false);
    else
        [~, gamma_breakthrough, ~] = computeGammaParameter3(effective.Swr, effective.Snwr, [], sw_front, ...
            fw_func, dfw_dsw_func, mob_t_func, 'Verbose', false, 'Process', effective.Process, 'Compute', 'immiscible', 'PlotIntegrand', false);
    end
    
    % Determine breakthrough PVI
    PVI_breakthrough = 1 / slope;
    
    % Compute effective initial total mobility
    mob_total_initial = mob_t_func(effective.SwInitial);

    % Initialize predictions
    fw_effective = zeros(size(PVI_analytical));
    dp_effective = fw_effective;
    
    for i = 1:length(PVI_analytical)
        PVI = PVI_analytical(i);
        % Before breakthrough
        if PVI <= PVI_breakthrough
            xD = PVI * slope;        % flood front location
            sw_out = effective.SwInitial;   % outlet saturation = initial saturation
            
            % Compute effective total mobility
            mob_total_effective_reciprocal = (xD/ gamma_breakthrough) + (1 / mob_total_initial) * (1 - xD);
            mob_total_effective = 1 ./ mob_total_effective_reciprocal;
            
        else
            % After breakthrough: solve f'(sw_out) = 1/PVI
            sw_out = fminbnd(@(s) abs(dfw_dsw_func(s) - 1/PVI), sw_min, sw_max);
            
            % Compute gamma parameter before breakthrough
            if effective.Miscible
                [mob_total_effective, ~, ~] = computeGammaParameter3(effective.Swr, effective.Snwr, sw_out, [], ...
                    fw_func, dfw_dsw_func, mob_t_func, 'Verbose', false, 'Process', effective.Process, 'Compute', 'miscible', 'PlotIntegrand', false);
            else
                [~, mob_total_effective, ~] = computeGammaParameter3(effective.Swr, effective.Snwr, [], sw_out, ...
                    fw_func, dfw_dsw_func, mob_t_func, 'Verbose', false, 'Process', effective.Process, 'Compute', 'immiscible', 'PlotIntegrand', false);
            end
            
        end
        
        % fractional flow 
        fw_effective(i) = fw_func(sw_out);
        % pressure drop
        % Compute flow capacity
        flow_capacity = effective.perm_abs * effective.thickness * mob_total_effective;
        % Compute pressure drop
        dp_effective(i) = 1 / flow_capacity;
    end
    
end

function Error = objective(theta, PVI_analytical, dp_analytical, fw_analytical, mu_w, mu_nw, effective)
% OBJECTIVE Sum of squared normalized errors in outlet fractional flow and
% pressure drop between the layered analytical solution and the
% equivalent homogeneous model (Eq. 7 of the manuscript).
    [dp_effective, fw_effective] = forwardHomogenized(theta, PVI_analytical, mu_w, mu_nw, effective);

    % Normalized errors (avoid division by zero)
    relative_error_fw = (fw_effective - fw_analytical) ./ max(fw_analytical);
    relative_error_dp = (dp_effective - dp_analytical) ./ max(dp_analytical);
    
    % Sum of squares
    Error = sum((relative_error_fw.^2) + (relative_error_dp.^2));

end