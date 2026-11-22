% =========================================================================
%% WORKFLOW FOR MULTILAYER SIMULATION
% =========================================================================

% DESCRIPTION:
%   Analytical solution for viscous-dominated two-phase displacement in
%   non-communicating stratified (layered) media, extended to semi-miscible
%   CO2-brine flow by combining the layered solution of Youssef and
%   Matthai (2022) with the semi-miscible fractional-flow theory of
%   Noh et al. (2007). Handles miscible (semi-miscible) and immiscible
%   cases and returns flood-front positions, outlet fractional flow,
%   pressure drop, free and dissolved volumes, and the dissolved-flux
%   correction factor alpha (miscible only).
%
%   The inputs below reproduce the Otway CRC-3 case of the accompanying
%   manuscript (Tables 2 and 3). Results are saved to results_file, which
%   is read by kr_inverse.m (upscaled relative permeability) and
%   soluble_flux.m.
%
% UNITS
%   Only ratios of thickness, permeability, viscosity and the two
%   mass/volume quantities enter the solution, so any consistent set of
%   units can be used. The Otway values are: thickness [cm],
%   permeability [mD], viscosity [cP], mass/volume [kg/m^3].

% INPUTS
% This code takes inputs as three patches: 1, 2 and 3. 
% Patch 1: geometrical and single phase parameters
% Patch 2: two phase parameters
% Patch 3: process and solution controls
% The user is required to change these patches to suit his case.
% The inputs shall be consistent between different patches not to get an
% error


close all;
clear;
clc;
warning off;
% =========================================================================
%% PARSE INPUTS
% =========================================================================

% --- INPUTS: Patch 1 (Geometrical and Single-Phase Parameters) ---
% SYNTAX:
%   layers - Structure array with fields:
%                      .porosity, .perm_abs, .thickness

% Layer 1: Very high permeability
layers(1).thickness = 31.1;      % Net flowing thickness [cm]
layers(1).porosity  = 0.214;     % Porosity [fraction]
layers(1).perm_abs  = 157.0;     % Absolute permeability [mD]

% Layer 2: High permeability
layers(2).thickness = 244.8;     % [cm]
layers(2).porosity  = 0.218;     % [-]
layers(2).perm_abs  = 129.3;     % [mD]

% Layer 3: Intermediate permeability
layers(3).thickness = 107.1;
layers(3).porosity  = 0.214;
layers(3).perm_abs  = 71.1;

% Layer 4: Low permeability
layers(4).thickness = 91.8;
layers(4).porosity  = 0.199;
layers(4).perm_abs  = 28.4;

% Layer 5: Very low permeability
layers(5).thickness = 8.0;
layers(5).porosity  = 0.195;
layers(5).perm_abs  = 8.3;
% =========================================================================

% --- INPUTS: Patch 2 (Two-Phase Parameters) ---

% SYNTAX:
%   mu_w      - Viscosity of wetting (aqueous) phase [any unit, same as mu_nw]
%   mu_nw     - Viscosity of non-wetting (CO2-rich) phase
%   models      - Cell array of strings specifying correlation models
%                       'Brooks-Corey' or 'BC' - Brooks-Corey model
%                       'Van-Genuchten' or 'VG' - Van Genuchten model
%                       'LET' - LET correlation model
%                       'Chierici' - Chierici model
%   params      - Cell array of structure containing parameters of selected
%                       model; e.g.,  params{1} = struct('lambda', 2.0);
%                       params{2} = struct('m', 0.5, 'n', 2.0); params{3} =
%                       struct('Lw',1.0,'Ew',1.5,'Tw',2.0,'Lnw',1.2,'Enw',1.0,'Tnw',1.8);
%                       Refere to "relativePermeability.m" for  more
%                       detials
%   Swr      - Residual water saturation (scalar or vector for multiple layers)
%   Snwr      - Residual non-wettong phase saturation (scalar or vector for multiple layers)
%   SwInitial      - Initial water saturation (scalar or vector for multiple layers)
%   Krw_end      - Endpoint relative permeability of wetting phase (scalar or vector for multiple layers)
%   Krnw_end      - Endpoint relative permeability of non-wetting phase (scalar or vector for multiple layers)

mu_w  = 0.482;                          % Viscosity of brine [cP]
mu_nw = 0.0409;                         % Viscosity of supercritical CO2 [cP]

models = {'BC', 'BC', 'BC', 'BC', 'BC'};  % Brooks-Corey model

Swr   = [0.2; 0.2; 0.2; 0.2; 0.2];     % Residual water saturation
Snwr  = [0.0; 0.0; 0.0; 0.0; 0.0];     % Residual non-wetting saturation
Krw_end  = [1.0; 1.0; 1.0; 1.0; 1.0];  % Water endpoint relative permeability
Krnw_end = [0.9; 0.9; 0.9; 0.9; 0.9];  % Non-wetting endpoint relative permeability
SwInitial = [1.0; 1.0; 1.0; 1.0; 1.0]; % Initial water saturation

% Model-specific parameters (e.g., lambda for Brooks-Corey), one cell per layer
params = cell(1, numel(layers));
for i = 1:numel(layers)
    params{i} = struct('lambda', 2.0);
end
% =========================================================================

% --- INPUTS: Patch 3 (Process and Simulation Controls) ---
% SYNTAX:
%   'Process'    -  String varaible determing if process is drainage or
%                        imbibition (Drainage, drainage, Imbibition or imbibition, default: imbibition)
%   solubility_into_displaced_phase  -  Mass of the displacing component (CO2) dissolved per unit
%                                       volume of the displaced (aqueous) phase, i.e. chi*rho_aq
%                                       [mass/volume]. Set to 0 for immiscible flow.
%   solubility_into_displacing_phase -  Mass of the displacing component per unit volume of the
%                                       displacing phase, i.e. the CO2-rich phase density rho_sc
%                                       when H2O vaporization is neglected [mass/volume].
%                                       The Welge tangent then starts at s = 1 + D, with
%                                       D = chi*rho_aq / (rho_sc - chi*rho_aq)  (Noh et al., 2007).
%   Resolution   -  Number of points used to discretize each stage
%   case_name    -  String used as a title of the case
%   results_file -  .mat file where results are saved (read by kr_inverse.m and soluble_flux.m)
%
%   NOTE: the miscible option is implemented for drainage (CO2 displacing brine,
%   SwInitial = 1). Miscible imbibition is not supported.

Process = 'drainage';                  % 'drainage' or 'imbibition'
solubility_into_displaced_phase  = 37.241;   % chi*rho_aq = 0.0365 x 1020.3 [kg/m^3] (Table 3)
solubility_into_displacing_phase = 545.7;    % rho_sc [kg/m^3] (Table 3)
Resolution   = 50;                           % Points per stage
case_name    = "Otway CRC-3 (semi-miscible)";
results_file = 'Otway_miscible_results.mat';

% =========================================================================
%% PRE-SIMULATION CALCULATIONS
% =========================================================================

% Validate layer properties and size consistency of per-layer inputs
n_layers = length(layers);
validateLayerProperties(layers, n_layers);
validateLayerInputs(n_layers, models, params, Swr, Snwr, Krw_end, Krnw_end, SwInitial);
Swr = Swr(:); Snwr = Snwr(:); Krw_end = Krw_end(:); Krnw_end = Krnw_end(:); SwInitial = SwInitial(:);

% Assign endpoint parameters to structure
for i = 1:n_layers
    params{i}.Swr = Swr(i);
    params{i}.Snwr = Snwr(i);
    params{i}.Krw_end = Krw_end(i);
    params{i}.Krnw_end = Krnw_end(i);
end

% Compute relative permeability and fractional flow functions
fprintf('\nComputing relative permeabilites and fractional flow functions...\n');
[krw_funcs, krnw_funcs, mob_w_funcs, mob_nw_funcs, mob_t_funcs, fw_funcs, dfw_dsw_funcs] = ...
    relativePermeability(models, params, mu_w, mu_nw);

% Compute total mobility and mobility ratios
mob_total_initial = zeros(n_layers, 1);
mobility_ratio    = zeros(n_layers, 1);
for i = 1:n_layers
    mob_total_initial(i) = mob_t_funcs{i}(SwInitial(i));
    if strcmpi(Process, 'drainage')
        mob_total_endpoint = mob_t_funcs{i}(Swr(i));
    else
        mob_total_endpoint = mob_t_funcs{i}(1 - Snwr(i));
    end
    mobility_ratio(i) = mob_total_endpoint / mob_total_initial(i);
end

% Determine miscibility state and tangent start points
Miscible = (solubility_into_displaced_phase > 0);
if Miscible && ~strcmpi(Process, 'drainage')
    error('MultiLayer:Unsupported', 'The miscible option is implemented for drainage only.');
end
if Miscible
    Sw_tanget_start = 1 - (solubility_into_displaced_phase / (solubility_into_displaced_phase - solubility_into_displacing_phase));
    Sw_tanget_start = Sw_tanget_start * ones(n_layers, 1);
else
    Sw_tanget_start = SwInitial;
end

% Compute shock front saturations
fprintf('\nComputing shock front saturations & gamma parameter...\n');
[sw_front, fw_front, slope] = findShockFront(Swr, Snwr, fw_funcs, dfw_dsw_funcs, ...
    Sw_tanget_start, 'Process', Process, 'Miscible', Miscible, 'Verbose', true, 'Plot', false);

% =========================================================================
%% FIRST STAGE SIMULATION
% =========================================================================

% Compute gamma parameter for the first stage
if Miscible
    [gamma, ~, ~] = computeGammaParameter3(Swr, Snwr, sw_front, [], ...
        fw_funcs, dfw_dsw_funcs, mob_t_funcs, 'Verbose', true, 'Process', Process, 'Compute', 'miscible');
else
    [~, gamma, ~] = computeGammaParameter3(Swr, Snwr, [], sw_front, ...
        fw_funcs, dfw_dsw_funcs, mob_t_funcs, 'Verbose', true, 'Process', Process, 'Compute', 'immiscible');
end

% Run first stage simulation
[stage1_results, layer_ordering] = simulateFirstStage(layers, gamma, sw_front, fw_front, ...
    slope, SwInitial, mob_total_initial, mobility_ratio, 'Process', Process, ...
    'Verbose', true, 'PlotResults', true, 'xD_resolution', Resolution, ...
    'case_name', case_name, 'solubility', solubility_into_displaced_phase);

fprintf('\nEnd of First Stage: === %s ', case_name);
if Miscible, fprintf('===Miscible flow===\n'); else, fprintf('===Immisble flow===\n'); end

% =========================================================================
%% SECOND STAGE SIMULATION
% =========================================================================

xD_end_first_stage = stage1_results.xD_all_layers(:, end);
solute_solvent_volume_ratio = stage1_results.solute_solvent_volume_ratio(end);

% Run second stage simulation (post-breakthrough)
[stage2_results, stage2_details] = simulateSecondStage3(layers, Swr, Snwr, gamma, sw_front, fw_front, ...
    slope, SwInitial, mob_total_initial, mobility_ratio, fw_funcs, dfw_dsw_funcs, mob_t_funcs, ...
    xD_end_first_stage, layer_ordering, solute_solvent_volume_ratio, 'Process', Process, ...
    'Verbose', true, 'n_points', Resolution, 'solubility', solubility_into_displaced_phase);

fprintf('\nEnd of Second Stage: === %s ', case_name);
if Miscible, fprintf('===Miscible flow===\n'); else, fprintf('===Immisble flow===\n'); end

% =========================================================================
%% COMBINED RESULTS
% =========================================================================
% This section merges the results from both stages of the simulation to 
% provide a continuous view of the storage performance and flow characteristics 
% from initial injection through post-breakthrough.

fprintf('\nProcessing combined results and generating plots...\n');

% 1. Aggregate numerical results into a single structure
combined.injected_pv   = [stage1_results.injected_pv, stage2_results.injected_pv];
combined.s_avg         = [stage1_results.avg_saturation, stage2_results.s_avg];
combined.free_volume   = [stage1_results.free_injected_volume, stage2_results.free_volume];
combined.soluble_volume = [stage1_results.soluble_injected_volume, stage2_results.soluble_volume];
combined.pressure_drop = [stage1_results.pressure_drop, stage2_results.pressure_drop];

% 2. Combine Fractional Flow (Fw) at the outlet
% Before breakthrough, the displaced phase fractional flow is 1 (for drainage)
if strcmpi(Process, 'drainage')
    fw_at_outlet_stage1 = ones(size(stage1_results.injected_pv));
else
    fw_at_outlet_stage1 = zeros(size(stage1_results.injected_pv));
end
combined.F_w = [fw_at_outlet_stage1, stage2_results.F_w];

% 3. Calculate Alpha Factor (Solubility efficiency)
% Alpha = average outlet solubility / solubility into the displaced phase
alpha_values_stage2 = stage2_results.solubility_outlet_average ./ (solubility_into_displaced_phase + eps);
combined.alpha_factor = [zeros(size(stage1_results.injected_pv)), alpha_values_stage2];

% Identify Breakthrough point for plot marking
pore_volume_at_breakthrough = stage1_results.injected_pv(end);

% --- FIGURE 1: Pore Volume Injected vs Storage Volumes ---
figure('Color', 'w', 'Name', 'Combined Storage Volumes', 'Position', [100, 100, 850, 600]);
hold on;
plot(combined.injected_pv, combined.free_volume * 100, 'b-', 'LineWidth', 3, 'DisplayName', 'Free Injected Phase');
plot(combined.injected_pv, combined.soluble_volume * 100, 'r-', 'LineWidth', 3, 'DisplayName', 'Soluble Injected Phase');
xline(pore_volume_at_breakthrough, '--k', 'Breakthrough', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Pore Volume Injected (PVI)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Storage Volume in Domain (%)', 'FontSize', 14, 'FontWeight', 'bold');
title('Combined Injected Phase Storage Volumes', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2);
grid on; box on;

% --- FIGURE 2: Pore Volume Injected vs Fw and Alpha Factor ---
figure('Color', 'w', 'Name', 'Combined Performance Characteristics', 'Position', [200, 100, 850, 600]);
yyaxis left
plot(combined.injected_pv, combined.F_w, 'b-', 'LineWidth', 3);
ylabel('Average Fractional Flow (F_w)', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 1.05]); set(gca, 'ycolor', 'b');

yyaxis right
plot(combined.injected_pv, combined.alpha_factor, 'r-', 'LineWidth', 3);
ylabel('Alpha Factor (\alpha)', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 1.05]); set(gca, 'ycolor', 'r');

xline(pore_volume_at_breakthrough, '--k', 'Breakthrough', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Pore Volume Injected (PVI)', 'FontSize', 14, 'FontWeight', 'bold');
title('Combined Performance Characteristics', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2);
grid on; box on;

fprintf('Combined analysis and plotting complete.\n');

% =========================================================================
%% SAVE RESULTS
% =========================================================================
% Variables required by kr_inverse.m, kr_inverse_BC.m and soluble_flux.m
save(results_file, 'combined', 'stage1_results', 'stage2_results', 'layers', ...
    'mu_w', 'mu_nw', 'models', 'params', 'Swr', 'Snwr', 'Krw_end', 'Krnw_end', ...
    'SwInitial', 'Process', 'Miscible', 'Sw_tanget_start', 'sw_front', 'fw_front', ...
    'solubility_into_displaced_phase', 'solubility_into_displacing_phase', ...
    'solute_solvent_volume_ratio', 'pore_volume_at_breakthrough');
fprintf('Results saved to %s\n', results_file);

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function validateLayerProperties(layers, n_layers)
% VALIDATELAYERPROPERTIES Check if layer structure has all required fields
    required_fields = {'thickness', 'porosity', 'perm_abs'};
    for i = 1:n_layers
        for j = 1:length(required_fields)
            if ~isfield(layers(i), required_fields{j})
                error('MultiLayer:MissingField', ...
                    'Layer %d missing field: %s', i, required_fields{j});
            end
        end
        
        % Validate values
        if layers(i).thickness <= 0
            error('MultiLayer:InvalidThickness', 'Layer %d: thickness <= 0', i);
        end
        if layers(i).porosity <= 0 || layers(i).porosity >= 1
            error('MultiLayer:InvalidPorosity', 'Layer %d: porosity must be (0,1)', i);
        end
        if layers(i).perm_abs <= 0
            error('MultiLayer:InvalidPermeability', 'Layer %d: perm_abs <= 0', i);
        end
    end
end

function validateLayerInputs(n_layers, models, params, varargin)
% VALIDATELAYERINPUTS Check that every per-layer input has one entry per layer
    if numel(models) ~= n_layers
        error('MultiLayer:SizeMismatch', ...
            '"models" has %d entries but there are %d layers', numel(models), n_layers);
    end
    if numel(params) ~= n_layers
        error('MultiLayer:SizeMismatch', ...
            '"params" has %d entries but there are %d layers', numel(params), n_layers);
    end
    names = {'Swr', 'Snwr', 'Krw_end', 'Krnw_end', 'SwInitial'};
    for k = 1:numel(varargin)
        if numel(varargin{k}) ~= n_layers
            error('MultiLayer:SizeMismatch', ...
                '"%s" has %d entries but there are %d layers', ...
                names{k}, numel(varargin{k}), n_layers);
        end
    end
end
