function [stage2_results, details] = simulateSecondStage3(layer_properties, Swr, Snwr, gamma, sw_front, fw_front,  ...
                                                  slope, SwInitial, mob_total_initial, mobility_ratio, fw_func, dfw_dsw_func, mob_t_func, ...
                                                  xD_end_first_stage, layer_ordering, solute_solvent_volume_ratio, varargin)
% SIMULATESECONDSTAGE Simulate second stage of unsteady state flooding for layered reservoir
%
% DESCRIPTION:
%   Computes flood front positions, outlet saturations, and reservoir performance
%   parameters for the second stage region (after breakthrough). This function builds on first stage
%   results and handles both pre-breakthrough and post-breakthrough regimes of each layer after breakthrough of fastest layer.
%
% INPUTS:
%   layer_properties - Structure array with fields:
%                      .porosity, .perm_abs, .thickness
%   Swr            - Residual water saturation [n_layers × 1]
%   Snwr           - Residual non-wetting saturation [n_layers × 1]
%   gamma      - Gamma parameter [n_layers × 1]
%   sw_front         - Shock front saturation for each layer [n_layers × 1]
%   fw_front         - Shock fractional flow for each layer [n_layers × 1]
%   slope      - dfw/dsw at front [n_layers × 1]
%   SwInitial            - Initial water saturation [n_layers × 1]
%   mob_total_initial   - Total mobility at initial saturation [n_layers × 1]
%   mob_ratio      - Endpoints mobility ratio
%   fw_func          - Cell array of fractional flow functions: fw(sw)
%   dfw_dsw_func     - Cell array of derivative functions: dfw_dsw(sw)
%   mob_t_func       - Cell array of total mobility functions: λ_T(sw)
%   xD_end_first_stage            - Flood front dimensionless distance at end of first
%                                               stage [n_layers × 1].
%   layer_ordering            - Ordering of layers from fastest to slowest
%                                          based on results obtained from
%                                          first stage [n_layers × 1].
%   solute_solvent_volume_ratio            - Dissolved volume relative to
%                                                               mixture volume after dissolution 
%
% OPTIONAL PARAMETERS:
%   'Process'        - 'drainage' or 'imbibition' (default: 'imbibition')
%   'solubility'        - Solubility for miscible cases (mass / volume).
%                                 (default: 0 assuming immiscible flow).
%   'n_points'       - Number of discretization points (default: 200)
%   'Tolerance'      - Tolerance for root finding (default: 1e-8)
%   'Verbose'        - Display progress (default: true)
%   'Plot'           - Plot results (default: true)
%   'MaxIter'        - Maximum iterations for bisection (default: 100)
%   'CombinePlots'   - Combine multiple plots (default: true)
%
% OUTPUTS:
%   results - Structure containing:
%             .sw_outlet_fastest - Outlet saturation in fastest layer
%             .xD_layers         - Dimensionless positions for all layers
%             .sw_outlet_layers  - Outlet saturations for all layers
%             .breakthrough_sw   - Breakthrough saturations for each layer
%             .s_avg             - Average saturation
%             .free_volume          - Free pore volume fraction of injected phase
%             .soluble_co2volume       - Soluble pore volume fraction of injected phase
%             .fastest_layer_idx - Index of fastest layer
%             .average_solubility          - Average solubility of model
%   details - Detailed intermediate calculations
%
% AUTHOR: AbdAllah AbdRabElnaby Youssef ElSaieed
% DATE: 04 February 2026
% =========================================================================
% PARSE INPUTS
% =========================================================================
% Default parameters
defaults.Process = 'imbibition';
defaults.solubility = 0;
defaults.n_points = 200;
defaults.Tolerance = 1e-8;
defaults.Verbose = true;
defaults.Plot = true;
defaults.MaxIter = 100;
defaults.CombinePlots = true;
% Parse name-value pairs
p = inputParser;
addRequired(p, 'layer_properties');
addRequired(p, 'Swr',  @(x) isnumeric(x));
addRequired(p, 'Snwr',  @(x) isnumeric(x));
addRequired(p, 'gamma', @(x) isnumeric(x) );
addRequired(p, 'sw_front', @(x) isnumeric(x));
addRequired(p, 'fw_front', @(x) isnumeric(x));
addRequired(p, 'slope', @(x) isnumeric(x));
addRequired(p, 'mob_total_initial', @(x) isnumeric(x));
addRequired(p, 'mobility_ratio', @(x) isnumeric(x));
addRequired(p, 'fw_func');
addRequired(p, 'dfw_dsw_func');
addRequired(p, 'mob_t_func');
addRequired(p, 'xD_end_first_stage', @(x) isnumeric(x));
addRequired(p, 'layer_ordering', @(x) isPositiveIntegerArray(x));
addRequired(p, 'solute_solvent_volume_ratio', @(x) isnumeric(x));
fields = fieldnames(defaults);
for i = 1:length(fields)
    addParameter(p, fields{i}, defaults.(fields{i}));
end
parse(p, layer_properties, Swr, Snwr, gamma, sw_front, fw_front, ...
          slope, mob_total_initial, mobility_ratio, fw_func, dfw_dsw_func, mob_t_func, ...
          xD_end_first_stage, layer_ordering, solute_solvent_volume_ratio, varargin{:});
params = p.Results;
% Validate required inputs
if min(Swr) < 0 || max(Swr) >= 1 || min(Snwr) < 0 || max(Snwr) >= 1 || max(Swr + Snwr) >= 1 || isempty(Swr) || isempty(Snwr)
    error('simulateSecondStage:InvalidResiduals', ...
           ['Input Swr & Snwr must be between zero and 1, ' ...
              '& their sum less than or equal 1']);
end
% Validate Process parameter
valid_process = {'drainage', 'imbibition', 'Drainage', 'Imbibition'};
params.Process = validatestring(params.Process, valid_process);
if strcmp(params.Process, 'Drainage') || strcmp(params.Process, 'drainage')
    params.Process = 'drainage';
else
    params.Process = 'imbibition';
end
% Validate flood front position
if any(xD_end_first_stage > 1 + 1e-12) || any(xD_end_first_stage <= 0)
    error('simulateSecondStage:InvalidFllodFrontPosition', ...
          'Flood front position at end of first stage must be more than 0 and not exceeding 1');   
    %%% ADD UNIQUE CHECK TO AVOID DUPLICATED LAYERS
end
% Validate total mobility at initial saturation
if any(mob_total_initial <= 0)
    error('simulateSecondStage:InvalidTotalMobilityAtInitialSaturation', ...
        'Total mobility at initial saturation must be more than 0');
end
% Check if solubility is valid
solubility = params.solubility;
if ~isnumeric(solubility)
        error('simulateSeondStage:UnvalidInput', ...
        'solubility must be scalar');
elseif solubility < 0     % immiscible flow
    solubility = 0;
end
% Check if solute-solvent volume ratio is valid
if ~isscalar(solute_solvent_volume_ratio) || solute_solvent_volume_ratio < 0
        error('simulateSeondStage:UnvalidInput', ...
        'solute-solvent volume ratio must be scalar more than or equal to zero');
end
% =========================================================================
% SETUP AND INITIALIZATION
% =========================================================================
n_layers = length(layer_properties);
% Extract layer properties
porosity = [layer_properties.porosity]';
thickness = [layer_properties.thickness]';
perm_abs = [layer_properties.perm_abs]';
% Compute first parameter: γ / (λ_T0 - γ)
first_param = gamma ./ (mob_total_initial - gamma);  % γ / (λ_T0 - γ)
% Compute second parameter 1: k * f' * λ_T0 / φ
second_param_1 = perm_abs .* mob_total_initial .* slope ./ porosity;      % k * f' * λ_T0 / φ
% Compute second parameter 2: k / φ
second_param_2 = perm_abs ./ porosity;      % k  / φ
% Determine fastest layer
if isempty(layer_ordering)
    error('simulateSecondStage:InvalidLayerOrdering', ...
              'Layer ordering must be non empty, positive integer array composed of %d elements', n_layers);
else
    layer_ordering = layer_ordering(:);     % transfering to column vector
    fastest_idx = layer_ordering(1);
end
% Ensure arrays are column vectors
Swr = Swr(:);
Snwr = Snwr(:);
sw_front = sw_front(:);
% Determine saturation bounds based on process
if strcmp(params.Process, 'drainage')
    sw_max = sw_front;  % For drainage, maximum is at front
    sw_min = Swr;
else
    sw_max = 1 - Snwr;  % For imbibition, maximum is 1-Snwr
    sw_min = sw_front;
end
if params.Verbose
    fprintf('\n=== Simulating Second Stage ===\n');
    fprintf('Number of layers: %d\n', n_layers);
    fprintf('Process: %s\n', params.Process);
    fprintf('Fastest layer: %d \n', fastest_idx);
    fprintf('Outlet saturation range of fastest layer: [%.4f, %.4f]\n', sw_min(fastest_idx), sw_max(fastest_idx));
    fprintf('\n');
end
% =========================================================================
% COMPUTE BREAKTHROUGH SATURATIONS FOR EACH LAYER
% =========================================================================
if params.Verbose
    fprintf('Computing breakthrough saturations...\n');
end
% Get alpha (double integral term) 
alpha = @(sw_star,layer_idx) computeDoubleIntegralTerm(Swr(layer_idx), Snwr(layer_idx), sw_front(layer_idx), sw_front(layer_idx), sw_star, ...
                              fw_func, dfw_dsw_func, mob_t_func, layer_idx, 'Process', params.Process, 'IntegrationMethod', 'trapz');   % Double integral parameter of fastest layer
% Breakthrough occurs when xD = 1 (front reaches outlet)
% Solve for sw_outlet_fastest where xD_layer_i = 1
breakthrough_sw = zeros(n_layers, 1);
breakthrough_converged = false(n_layers, 1);
second_param_2_fastest_layer = second_param_2(fastest_idx);
parfor layer = 1:n_layers
    if layer == fastest_idx
        breakthrough_sw(layer) = sw_front(fastest_idx);
        breakthrough_converged(layer) = true;
        continue;
    end
    
    if params.Verbose
        fprintf('  Layer %d ', layer);
    end
    
    xD_end_first_stage_layer = xD_end_first_stage(layer);
    % Define function to find root of: xD(sw_outlet_fastest) 
    f_breakthrough = @(sw_fastest) computeFrontPosition(alpha(sw_fastest, fastest_idx), first_param(layer), second_param_1(layer), ...
                                                       second_param_2_fastest_layer, mobility_ratio(layer));
    f_breakthrough = @(sw_fastest) xD_end_first_stage_layer + f_breakthrough(sw_fastest);
    f_breakthrough_target = @(sw_fastest) abs(f_breakthrough(sw_fastest) - 1);
    % Set search bounds
    lower_bound = sw_min(fastest_idx);
    upper_bound = sw_max(fastest_idx);
    
    % Use fminbnd to find breakthrough saturation
    [sw_bt, fval_bt, exit_flag] = fminbnd(f_breakthrough_target, lower_bound , upper_bound);

    % fminbnd is a local minimizer; when f_breakthrough is infinite over part
    % of the interval it can stop on that plateau. If the minimum found does
    % not satisfy xD = 1, locate the root by grid bracketing + fzero instead.
    if exit_flag <= 0 || ~(fval_bt < 1e-3)
        [sw_bt_alt, ok_alt] = bracketedRoot(@(s) f_breakthrough(s) - 1, lower_bound, upper_bound);
        if ok_alt
            sw_bt = sw_bt_alt;
            exit_flag = 1;
        end
    end
    
    if exit_flag > 0
        breakthrough_sw(layer) = sw_bt;
        breakthrough_converged(layer) = true;
        
        if params.Verbose
            fprintf('sw of fastest layer @ bt of Layer %d = %.4f\n', layer, sw_bt);
        end
    else
        warning('simulateSecondStage:BreakthroughNotConverged', ...
                'Breakthrough not found for layer %d', layer);
        breakthrough_sw(layer) = upper_bound;
        breakthrough_converged(layer) = false;
        
        if params.Verbose
            fprintf('failed to converge\n');
        end
    end
end
% =========================================================================
% CREATE SATURATION GRID FOR FASTEST LAYER
% =========================================================================
% Include breakthrough points in the grid
unique_bt_points = unique([sw_min(fastest_idx); breakthrough_sw; sw_max(fastest_idx)]);
if strcmp(params.Process, 'drainage')
    unique_bt_points = sort(unique_bt_points, 'descend');
else
    unique_bt_points = sort(unique_bt_points, 'ascend');
end
% Create adaptive grid
sw_segments = cell(length(unique_bt_points)-1, 1);
for i = 1:length(unique_bt_points)-1
    segment_start = unique_bt_points(i);
    segment_end = unique_bt_points(i+1);
    
    % Determine number of points based on segment length
    n_segment = max(3, round(params.n_points * abs(segment_end - segment_start) / ...
                            abs(sw_max(fastest_idx) - sw_min(fastest_idx))));
    
    sw_segments{i} = linspace(segment_start, segment_end, n_segment)';
end
% Combine segments keeping the duplicates at breakthrough points to capture steps
sw_outlet_fastest = vertcat(sw_segments{:});
% Identify second occurrences (post-BT states) for breakthrough points
is_second_occurrence = false(size(sw_outlet_fastest));
for counter = 2:length(sw_outlet_fastest)
    if sw_outlet_fastest(counter) == sw_outlet_fastest(counter-1)
        is_second_occurrence(counter) = true;
    end
end
n_points_total = length(sw_outlet_fastest);
if params.Verbose
    fprintf('Created saturation grid for fastest layer\n');
    fprintf('  Number of points: %d\n', n_points_total);
    fprintf('  Range: [%.4f, %.4f]\n', min(sw_outlet_fastest), max(sw_outlet_fastest));
    fprintf('\n');
end
% =========================================================================
% COMPUTE FRONT POSITIONS AND OUTLET SATURATIONS
% =========================================================================
if params.Verbose
    fprintf('Computing front positions and outlet saturations...\n');
end
% Initialize arrays
xD_layers = zeros(n_layers, n_points_total);      % Dimensionless positions
sw_outlet_layers = zeros(n_layers, n_points_total); % Outlet saturations
mob_total_effective = zeros(n_layers, n_points_total); % Effective total mobility
layer_status = zeros(n_layers, n_points_total);   % 0=pre-BT, 1=post-BT
% Get alpha of first layer at time of breakthrough of other layers.
alpha_bt = zeros(n_layers,1);
parfor layer = 1:n_layers
    if layer ~= fastest_idx
        alpha_bt(layer) = alpha(breakthrough_sw(layer), fastest_idx);
    end
end
% For each point in fastest layer saturation
parfor point_idx = 1:n_points_total
    sw_fastest = sw_outlet_fastest(point_idx);
    is_post = is_second_occurrence(point_idx);
    
    if params.Verbose && mod(point_idx, 50) == 0
        fprintf('  Processing point %d/%d: sw_fastest = %.4f\n', ...
                point_idx, n_points_total, sw_fastest);
    end
    
    alpha_fastest = alpha(sw_fastest, fastest_idx);
    % For each layer
    for layer = 1:n_layers
        if layer == fastest_idx
            % Fastest layer is always at outlet saturation
            % xD_layers(layer, point_idx) = 1;
            sw_outlet_layers(layer, point_idx) = sw_fastest;
            layer_status(layer, point_idx) = 1;  % Always post-BT
            % Compute effective total mobility
            if strcmp(params.Process, 'drainage') && sw_fastest == Swr(fastest_idx)
                mob_total_effective(layer,point_idx) = mob_t_func{layer}(Swr(fastest_idx));
            elseif strcmp(params.Process, 'imbibition') && sw_fastest == 1 - Snwr(fastest_idx)
                mob_total_effective(layer,point_idx) = mob_t_func{layer}(1 - Snwr(fastest_idx));
            elseif solubility == 0
                [~, mob_total_effective(layer,point_idx), ~]  = computeGammaParameter3(Swr(layer), Snwr(layer), [], sw_fastest, ...
                          fw_func{layer}, dfw_dsw_func{layer}, mob_t_func{layer}, 'Verbose', false, 'Process', params.Process, 'Compute' ,'immiscible', 'PlotIntegrand', false);
            else
                [mob_total_effective(layer,point_idx), ~, ~]  = computeGammaParameter3(Swr(layer), Snwr(layer), sw_fastest, [], ...
                          fw_func{layer}, dfw_dsw_func{layer}, mob_t_func{layer}, 'Verbose', false, 'Process', params.Process, 'Compute' ,'miscible', 'PlotIntegrand', false);
            end
            continue;
        end
        
        % Determine if layer has broken through at this sw_fastest
         if strcmp(params.Process, 'drainage')
            has_broken_through = (sw_fastest < breakthrough_sw(layer));
            if sw_fastest == Swr(fastest_idx)
                sw_outlet_layers(layer, point_idx) = Swr(layer);
                xD_layers(layer, point_idx) = 1 ;   % - xD_end_first_stage(layer);
                % Compute effective total mobility
                mob_total_effective(layer,point_idx) = mob_t_func{layer}(Swr(layer));
                layer_status(layer, point_idx) = 1;
                continue;
            end
         else
             has_broken_through = (sw_fastest > breakthrough_sw(layer));
             if sw_fastest == 1 - Snwr(fastest_idx)
                 sw_outlet_layers(layer, point_idx) = 1 - Snwr(layer);
                 xD_layers(layer, point_idx) = 1;    %  - xD_end_first_stage(layer);
                % Compute effective total mobility
                mob_total_effective(layer,point_idx) = mob_t_func{layer}(1 - Snwr(layer));
                 layer_status(layer, point_idx) = 1;
                 continue;
             end
         end

         % Treating jump at breakthrough
         if sw_fastest == breakthrough_sw(layer)
             if is_post
                 % Post-breakthrough state
                 sw_outlet_layers(layer, point_idx) = sw_front(layer);
                 layer_status(layer, point_idx) = 1;
             else
                 % Pre-breakthrough state at breakthrough point
                 sw_outlet_layers(layer, point_idx) = SwInitial(layer);                
                 layer_status(layer, point_idx) = 0;
             end
             xD_layers(layer, point_idx) = 1 ;
             % Compute effective total mobility at BT
             mob_total_effective(layer,point_idx) =  gamma(layer);
             continue;
         end
        
         first_param_layer = first_param(layer);
         second_param_1_layer = second_param_1(layer);
         mobility_ratio_layer = mobility_ratio(layer);
        if ~has_broken_through
            % Pre-breakthrough: compute front position
            xD_layers(layer, point_idx) = computeFrontPosition(alpha_fastest, first_param_layer, second_param_1_layer, ...
                                                       second_param_2_fastest_layer, mobility_ratio_layer);
            xD_layers(layer, point_idx) = xD_layers(layer, point_idx) + xD_end_first_stage(layer); 
            layer_status(layer, point_idx) = 0;  % Pre-breakthrough
            
            % Compute effective total mobility
            mob_total_effective_reciprocal = (xD_layers(layer, point_idx) / gamma(layer)) + (1 / mob_total_initial(layer)) * (1 - xD_layers(layer, point_idx));
            mob_total_effective(layer,point_idx) =  1 / mob_total_effective_reciprocal;
            sw_outlet_layers(layer, point_idx) = SwInitial(layer);
        else
            % Post-breakthrough: compute outlet saturation
            % Solve for sw_outlet 
            
            alpha_fastest_bt = alpha(sw_fastest, fastest_idx) - alpha_bt(layer);
            f_outlet_layer = @(sw_out) abs( alpha(sw_out, layer) - (second_param_2(layer) / second_param_2_fastest_layer) * alpha_fastest_bt);
            
            % Search bounds for outlet saturation
            if strcmp(params.Process, 'imbibition')
                lower_sw = sw_front(layer);
                upper_sw = sw_max(layer);
            else
                lower_sw = sw_min(layer);
                upper_sw = sw_front(layer);
            end
            
            % Use fminbnd method to find outlet saturation
            [sw_outlet, ~, exit_flag] = fminbnd(f_outlet_layer, lower_sw, upper_sw);
            
            if exit_flag > 0
                sw_outlet_layers(layer, point_idx) = sw_outlet;
                layer_status(layer, point_idx) = 1;  % Post-breakthrough
            else
                warning('simulateSecondStage:OutletSaturationNotConverged', ...
                        'Layer %d: Outlet saturation not found at sw_fastest = %.4f', ...
                        layer, sw_fastest);
                sw_outlet_layers(layer, point_idx) = upper_sw;  % Use upper bound
                layer_status(layer, point_idx) = 1;
            end
            % Compute effective total mobility
            if solubility == 0
                [~, mob_total_effective(layer,point_idx), ~]  = computeGammaParameter3(Swr(layer), Snwr(layer), [], sw_outlet_layers(layer, point_idx), ...
                          fw_func{layer}, dfw_dsw_func{layer}, mob_t_func{layer}, 'Verbose', false, 'Process', params.Process, 'Compute' ,'immiscible', 'PlotIntegrand', false);
            else
                [mob_total_effective(layer,point_idx), ~, ~]  = computeGammaParameter3(Swr(layer), Snwr(layer), sw_outlet_layers(layer, point_idx), [], ...
                          fw_func{layer}, dfw_dsw_func{layer}, mob_t_func{layer}, 'Verbose', false, 'Process', params.Process, 'Compute' ,'miscible', 'PlotIntegrand', false);
            end
            
        end
    end
end
xD_layers(layer_status==1) = 1;
% Ensure physical bounds
% xD_layers = max(0, min(xD_layers, 1));

% =========================================================================
% COMPUTE FLOW RATE PARTITIONING & PRESSURE DROP
% =========================================================================

% Compute layer flow capacities
flow_capacity = (perm_abs .* thickness) .* mob_total_effective;
% Compute flow rate ratio relative to fastest layer
layer_rate_realtive_to_fastest_layer = flow_capacity ./ flow_capacity(fastest_idx, :);
% Compute flow rate of fastest layer
fasest_layer_rate = 1 ./sum(layer_rate_realtive_to_fastest_layer);
% Compute flow rate in all layers
rate_layers = fasest_layer_rate .* layer_rate_realtive_to_fastest_layer;
% Compute pressure drop
pressure_drop = fasest_layer_rate ./ flow_capacity(fastest_idx, :);
% =========================================================================
% COMPUTE OUTLET SATURATION AND FRACTIONAL FLOW
% =========================================================================
if params.Verbose
    fprintf('\n--- Computing Average Saturations and Volumes of Second Stage ---\n');
end
% Compute average paramters at outlet
[sw_outlet_average, free_volume_outlet, fw_outlet_average, fw_outlet, dfw_dsw_outlet, solubility_outlet_average] = computeSliceAverage(sw_outlet_layers, rate_layers, layer_status, porosity , thickness, fw_func, dfw_dsw_func, solubility, params.Process);
% Compute domian average parameters 
[avg_saturation, free_volume, soluble_volume, total_storage_volume, injected_pv, average_solubility] = ...
    computeVolumes( xD_layers, sw_outlet_layers, fw_outlet, sw_front, fw_front, slope, dfw_dsw_outlet, layer_status, ...
   porosity, thickness, Swr, Snwr, SwInitial, solubility, solute_solvent_volume_ratio, params.Process);
% =========================================================================
% ORGANIZE RESULTS
% =========================================================================
stage2_results = struct();
stage2_results.sw_outlet_fastest = sw_outlet_fastest;
stage2_results.sw_outlet_average = sw_outlet_average;
stage2_results.solubility_outlet_average = solubility_outlet_average;
stage2_results.xD_layers = xD_layers;
stage2_results.sw_outlet_layers = sw_outlet_layers;
stage2_results.layer_status = layer_status;
stage2_results.breakthrough_sw = breakthrough_sw;
stage2_results.breakthrough_converged = breakthrough_converged;
stage2_results.s_avg = avg_saturation;
stage2_results.average_solubility = average_solubility;
stage2_results.free_volume = free_volume;
stage2_results.free_volume_outlet = free_volume_outlet;
stage2_results.soluble_volume = soluble_volume;
stage2_results.total_storage_volume = total_storage_volume;
stage2_results.injected_pv = injected_pv;
stage2_results.F_w = fw_outlet_average;
stage2_results.fastest_layer_idx = fastest_idx;
stage2_results.process = params.Process;
stage2_results.pressure_drop = pressure_drop;
stage2_results.rate_layers = rate_layers;
details = struct();
details.params = params;
details.layer_properties = layer_properties;
details.sw_front = sw_front;
details.gamma = gamma;
% =========================================================================
% VISUALIZATION
% =========================================================================
if params.Plot
    plotSecondStageResults(stage2_results, details, params);
end
if params.Verbose
    fprintf('\n=== Second Stage Simulation Complete ===\n');
    fprintf('Summary:\n');
    fprintf('  Fastest layer breakthrough: end sw of fastest layer = %.4f\n', breakthrough_sw(fastest_idx));
    for layer = 1:n_layers
        fprintf('  Layer %d breakthrough: end sw of fastest layer = %.4f\n', layer, breakthrough_sw(layer));
    end
    fprintf('\n');
end
end % Main function
% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function xD_layer_i = computeFrontPosition(alpha_fastest, first_param_i, second_param_1_i, ...
                                                       second_param_2_fastest_layer, mob_ratio_i)
% COMPUTELAYERPOSITION Compute dimensionless position for layer i
    
    % Extract parameters
    first_param_i = first_param_i(:);  % Ensure column vector [n_layers × 1]
    second_param_1_i = second_param_1_i(:);   % Ensure column vector [n_layers × 1]
    %%% NEED TO ADD OPERATOR THAT HANDLES A FUNCTION INSTEAD OF TRANSPOSE 
    alpha_fastest = alpha_fastest(:)';   
    
    % Compute using the equation from AbdAllah & Matthai (2022)
    term2 = (2 ./  first_param_i).* second_param_1_i / second_param_2_fastest_layer;
    
     radicand = 1 -  alpha_fastest .* term2;
    
    % Check for negative radicands
    negative_mask = radicand < 0;
    if any(negative_mask(:))
        warning('simulateSecondStage:NegativeRadicand', ...
                'Negative radicand detected in %d combinations. Using absolute value.', ...
                sum(negative_mask(:)));
        radicand(negative_mask) = inf; % abs(radicand(negative_mask));
    end
        
    % Choosing correct formula of position in layer i based on endpoint mobility ratio
    if all(mob_ratio_i ~= 1)
        xD_layer_i = first_param_i .* (-1 + sqrt(radicand));     
    elseif all(mob_ratio_i == 1)
        term2 = second_param_1_i / second_param_2_fastest_layer;
        xD_layer_i = alpha_fastest .* term2 ;
    else
        n_layers = length(first_param_i);
        n_xD_fastest = length(xD_fastest);
        xD_layer_i = zeros(n_layers , n_xD_fastest);
        for i = 1:n_layers
            if mob_ratio_i(i) ~= 1
                xD_layer_i(i,:) = first_param_i(i) .* (-1 + sqrt(radicand(i,:)));   
            else
                term2(i) = second_param_1_i(i) / second_param_2_fastest_layer;
                xD_layer_i(i,:) = alpha_fastest .* term2(i);
            end
        end
    end
    
    % Checking if there are some vlaues of flood front position are
    % non-physical (negative)
    % if any(xD_layer_i < 0)
    %     [non_physical_layer_idx, ~] = find(xD_layer_i < 0);
    %     non_physical_layer_idx = unique(non_physical_layer_idx);
    %     error('simulateSecondStage:NonPhysicalValues', ...
    %                   'Layer %d gives non physical flood front position', non_physical_layer_idx);
    % end
end
function [avg_sat_total, free_volume, soluble_volume, total_storage_volume, injected_pv, avg_solubility_total] = computeVolumes(...
    xD_all_layers, sw_outlet_all_layers, fw_outlet_all_layers, sw_front, fw_front, slope_front, slope_outlet, layer_status, porosity, thickness, Swr, Snwr, ...
    sw_initial, solubility, solute_solvent_volume_ratio, Process)
% COMPUTEVolumes Compute average saturations and volumes dissoloved and
% free
    
% Compute pore volumes 
    pore_volume = porosity .* thickness;     % normalized pore volume of each layer
    total_pore_volume = sum(pore_volume);  % total pore volume
    
% Compute injected volumes for layers that did not reach breakthrough
    injected_volume_stage_1 = (pore_volume ./ slope_front) .* xD_all_layers;     % This estimation is correct for immiscible and approximation for immiscible 
    injected_volume_stage_1(layer_status == 1) = 0;
% Compute injected volumes for layers that reached breakthrough
    injected_volume_stage_2 = pore_volume ./ (slope_outlet + 1e-10);     % This estimation is correct for immiscible and approximation for immiscible 
    injected_volume_stage_2(layer_status == 0) = 0;
% Compute combined injected volumes
    injected_volume = injected_volume_stage_1 + injected_volume_stage_2;
    total_injecetd_volume = sum(injected_volume);       
% Compute pore volume injected
    injected_pv = total_injecetd_volume / total_pore_volume;
    injected_pv(isnan(injected_pv)) = inf;

% Compute average saturation that did not reach breakthrough and correct
% injected pore volime at extreme point
    if strcmp(Process, 'drainage')
        fw_inlet = 0;
    else
        fw_inlet = 1;
    end
    sw_avg_behind = sw_front - (1 ./ slope_front) .* (fw_front - fw_inlet);     % average saturation behind front
    sw_avg_stage_1 = sw_avg_behind .* xD_all_layers + sw_initial .* (1 - xD_all_layers);   % average saturation in each layer
    sw_avg_stage_1(layer_status == 1) = 0;
    
% Compute average saturation that reached breakthrough
    sw_avg_stage_2 = sw_outlet_all_layers - (fw_outlet_all_layers - fw_inlet) ./ (slope_outlet + 1e-10);     % average saturation in each layer
    
    % Robust Safety Clipping for sw_avg_stage_2
    if strcmp(Process, 'drainage')
        % For drainage, Swr <= sw_avg <= sw_outlet
        sw_avg_stage_2 = min(max(sw_avg_stage_2, Swr), sw_outlet_all_layers);
        if injected_pv(end) > 1e10 || isinf(injected_pv(end))
            sw_avg_stage_2(:, end) = Swr;
        end
    else
        % For imbibition, sw_outlet <= sw_avg <= 1-Snwr
        sw_avg_stage_2 = max(min(sw_avg_stage_2, 1 - Snwr), sw_outlet_all_layers);
        if injected_pv(end) > 1e10 || isinf(injected_pv(end))
            sw_avg_stage_2(:, end) = 1 - Snwr;
        end
    end
    sw_avg_stage_2(layer_status == 0) = 0;
% Compute combined average saturation
    sw_avg = sw_avg_stage_1 + sw_avg_stage_2;
% Compute water volume
    water_volume = pore_volume .* sw_avg;
    water_volume_total = sum(water_volume);
% Compute model avearge saturation
    avg_sat_total = water_volume_total / total_pore_volume;

% Compute free phase and soluble volume  
   if strcmp(Process, 'drainage')
       free_volume = 1 - avg_sat_total;
        
       % Compute solvent volume for layers that did not reach breakthrough
       soluble_volume_stage_1 = sw_avg_behind .* xD_all_layers;
       soluble_volume_stage_1(layer_status == 1) = 0;
       % Compute solvent volume for layers that reached breakthrough
       soluble_volume_stage_2 = sw_avg_stage_2;
       soluble_volume_stage_2(layer_status == 0) = 0;
       
   else
       free_volume = avg_sat_total;
       % Compute solvent volume for layers that did not reach breakthrough
       soluble_volume_stage_1 = (1 - sw_avg_behind) .* xD_all_layers;
       soluble_volume_stage_1(layer_status == 1) = 0;
       % Compute solvent volume for layers that reached breakthrough
       soluble_volume_stage_2 = 1 - sw_avg_stage_2;
       soluble_volume_stage_2(layer_status == 0) = 0;
   end
   
   % Compute combined solvent volume
    solvent_volume = pore_volume .* (soluble_volume_stage_1 + soluble_volume_stage_2);
    
    % Compute soluble volume
    soluble_volume = sum(solute_solvent_volume_ratio * solvent_volume) / total_pore_volume;
    % Compute total storage volume
    total_storage_volume = free_volume + soluble_volume;
    % Compute soluble mass
    soluble_mass = sum(solubility * solvent_volume);
% Compute model average solubility
if solubility > 0
    if strcmp(Process, 'drainage')
        avg_solubility_total = soluble_mass ./ water_volume_total;
    else
        avg_solubility_total = soluble_mass ./ (total_pore_volume - water_volume_total);
    end
else
    avg_solubility_total = 0 * injected_pv;
end
end

function [avg_sat_slice, free_phase, avg_fw_slice, fw_slice,dfw_dsw_slice,  avg_solubility] = computeSliceAverage( ...
    sw_all_layers, rate_all_layers, layer_status, porosity , thickness, fw_func, dfw_dsw_func, solubility, Process)
% COMPUTESliceAverage Compute average saturations, volumes (dissolved
% and free), fractional flow, fractional flow derivative distribution and average solubility at certain slice
% Compute pore volumes 
    pore_volume = porosity .* thickness;     % normalized pore volume of each layer
    total_pore_volume = sum(pore_volume);  % total pore volume
    
% Compute water volume
     water_volume = pore_volume .* sw_all_layers;
     water_volume_total = sum(water_volume);
% Compute average saturation
    avg_sat_slice = water_volume_total / total_pore_volume;
% Compute fractional flow and fractional flow derivative
    [n_layers , n_points] = size(sw_all_layers);
    fw_slice = zeros(n_layers, n_points);
    dfw_dsw_slice = zeros(n_layers, n_points);
    for layer = 1:n_layers
        fw_slice(layer, :) = fw_func{layer}(sw_all_layers(layer, :));
        dfw_dsw_slice(layer, :) = dfw_dsw_func{layer}(sw_all_layers(layer, :));
    end
% Compute water flow rate
    qw_slice = fw_slice .* rate_all_layers;
    avg_fw_slice  = sum(qw_slice);
% Compute free phase volume
   if strcmp(Process, 'drainage')
       free_phase = 1 - avg_sat_slice;
   else
       free_phase = avg_sat_slice;
   end
% Compute average solubility
% Solubility is only calcualted for layers that brokethrough
    if solubility > 0
        if strcmp(Process, 'drainage')
             soluble_phase_rate_slice = layer_status .* qw_slice;
             avg_solubility =solubility * sum(soluble_phase_rate_slice) ./ avg_fw_slice;
        else
            qnw_slice = (1 - fw_slice) .* rate_all_layers;
            soluble_phase_rate_slice = layer_status .* qnw_slice;
            avg_solubility = solubility * sum(soluble_phase_rate_slice) ./ (1 - avg_fw_slice);
        end
    else
        avg_solubility = 0 * avg_sat_slice;
    end
end
function plotSecondStageResults(results, details, params)
% PLOTSECONDSTAGERESULTS Create comprehensive visualization
    
    n_layers = length(details.sw_front);
    fastest_idx = results.fastest_layer_idx;
    
    % Determine layout
    if params.CombinePlots && n_layers <= 6
        % Combined figure
        createCombinedSecondStageFigure(results, details, params);
    else
        % Separate figures
        createSeparateSecondStageFigures(results, details, params);
    end
end
function createCombinedSecondStageFigure(results, details, params)
% CREATECOMBINEDSECONDSTAGEFIGURE Combined visualization
    
    n_layers = length(details.sw_front);
    
    fig = figure('Position', [100, 100, 1800, 1200], 'Color', 'white');
    
    % Subplot 1: Front positions vs outlet saturation
    subplot(2, 3, 1);
    hold on;
    
    colors = lines(n_layers);
    
    for layer = 1:n_layers
        if layer == results.fastest_layer_idx
            style = '--';
            linewidth = 3;
        else
            style = '-';
            linewidth = 1.5;
        end
        
        plot(results.sw_outlet_fastest, results.xD_layers(layer, :), ...
             style, 'Color', colors(layer, :), 'LineWidth', linewidth);
        
        % Mark breakthrough point
        if results.breakthrough_converged(layer) && layer ~= results.fastest_layer_idx
            x_bt = results.breakthrough_sw(layer);
            y_bt = 1;  % At breakthrough
            plot(x_bt, y_bt, 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', colors(layer, :), ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Dimensionless Position (x_D)');
    title('Front Positions vs Fastest Layer Saturation');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    box on;
    
    % Subplot 2: Outlet saturations
    subplot(2, 3, 2);
    hold on;
    
    for layer = 1:n_layers
        if layer == results.fastest_layer_idx
            style = '--';
            linewidth = 3;
        else
            style = '-';
            linewidth = 1.5;
        end
        
        plot(results.sw_outlet_fastest, results.sw_outlet_layers(layer, :), ...
             style, 'Color', colors(layer, :), 'LineWidth', linewidth);
    end
    
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Outlet Saturation in Each Layer (s_w)');
    title('Layer Outlet Saturations');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    box on;
    
    % Subplot 3: Average saturation
    subplot(2, 3, 3);
    plot(results.sw_outlet_fastest, results.s_avg, 'b-', 'LineWidth', 2);
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Average Saturation (s_{avg})');
    title('Average Reservoir Saturation');
    grid on;
    box on;
    
    % % Subplot 4: Recovery factor
    % subplot(2, 3, 4);
    % plot(results.sw_outlet_fastest, results.recovery * 100, 'r-', 'LineWidth', 2);
    % xlabel('Outlet Saturation in Fastest Layer (s_w)');
    % ylabel('Recovery Factor (%)');
    % title('Recovery Factor');
    % grid on;
    % box on;
    
    % Subplot 5: Free CO2
    subplot(2, 3, 5);
    plot(results.sw_outlet_fastest, results.free_volume * 100, 'g-', 'LineWidth', 2);
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Injected phase free pore volume (%)');
    title('Free pore volume Fraction');
    grid on;
    box on;
    
    % Subplot 6: Summary
    subplot(2, 3, 6);
    axis off;
    
    % summary_text = sprintf(['Second Stage Results\n\n' ...
    %                        'Process: %s\n' ...
    %                        'Fastest Layer: %s\n' ...
    %                        'Number of Layers: %d\n\n' ...
    %                        'Breakthrough Order:\n'], ...
    %                       results.process, ...
    %                       results.layer_names{results.fastest_layer_idx}, ...
    %                       n_layers);
    % 
    % for i = 1:min(5, length(results.breakthrough_sw))
    %     [~, order] = sort(results.breakthrough_sw);
    %     layer = order(i);
    %     if results.breakthrough_converged(layer)
    %         summary_text = sprintf('%s%d. %s: s_w = %.4f\n', ...
    %                               summary_text, i, ...
    %                               results.layer_names{layer}, ...
    %                               results.breakthrough_sw(layer));
    %     end
    % end
    % 
    % summary_text = sprintf('%s\nFinal Values:\n', summary_text);
    % summary_text = sprintf('%s  Avg. Saturation: %.4f\n', ...
    %                       summary_text, results.s_avg(end));
    % summary_text = sprintf('%s  Recovery: %.2f%%\n', ...
    %                       summary_text, results.recovery(end) * 100);
    % summary_text = sprintf('%s  Free CO2: %.2f%%', ...
    %                       summary_text, results.free_co2(end) * 100);
    % 
    % text(0.1, 0.5, summary_text, 'FontSize', 10, ...
    %      'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    % 
    % sgtitle('Second Stage Analysis - Combined Results', ...
    %         'FontSize', 16, 'FontWeight', 'bold');
end
function createSeparateSecondStageFigures(results, details, params)
% CREATESEPARATESECONDSTAGEFIGURES Separate visualizations
    
    % Figure 1: Front positions
    fig1 = figure('Position', [100, 100, 1200, 800], 'Color', 'white');
    hold on;
    
    colors = lines(length(results.layer_names));
    
    for layer = 1:length(results.layer_names)
        if layer == results.fastest_layer_idx
            style = '--';
            linewidth = 3;
        else
            style = '-';
            linewidth = 1.5;
        end
        
        plot(results.sw_outlet_fastest, results.xD_layers(layer, :), ...
             style, 'Color', colors(layer, :), 'LineWidth', linewidth);
    end
    
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Dimensionless Position (x_D)');
    title('Front Positions vs Fastest Layer Saturation');
    legend('Location', 'best');
    grid on;
    box on;
    
    % Figure 2: Performance parameters
    fig2 = figure('Position', [200, 100, 1400, 600], 'Color', 'white');
    
    subplot(1, 3, 1);
    plot(results.sw_outlet_fastest, results.s_avg, 'b-', 'LineWidth', 2);
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Average Saturation (s_{avg})');
    title('Average Saturation');
    grid on;
    box on;
    
    subplot(1, 3, 2);
    plot(results.sw_outlet_fastest, results.recovery * 100, 'r-', 'LineWidth', 2);
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Recovery Factor (%)');
    title('Recovery Factor');
    grid on;
    box on;
    
    subplot(1, 3, 3);
    plot(results.sw_outlet_fastest, results.free_co2 * 100, 'g-', 'LineWidth', 2);
    xlabel('Outlet Saturation in Fastest Layer (s_w)');
    ylabel('Free CO_2 (%)');
    title('Free CO_2 Fraction');
    grid on;
    box on;
    
    sgtitle('Second Stage Performance Parameters', ...
            'FontSize', 14, 'FontWeight', 'bold');
end


function [root, ok] = bracketedRoot(g, lb, ub)
% BRACKETEDROOT Find a root of g on [lb, ub]: evaluate g on a grid, take the
% first sign change between finite real values, and refine it with fzero.
    n  = 201;
    s  = linspace(lb, ub, n);
    gv = nan(1, n);
    for k = 1:n
        v = g(s(k));
        if isreal(v) && isfinite(v)
            gv(k) = v;
        end
    end
    root = NaN; ok = false;
    for k = 1:n-1
        if ~isnan(gv(k)) && ~isnan(gv(k+1)) && sign(gv(k)) ~= sign(gv(k+1))
            root = fzero(g, [s(k), s(k+1)]);
            ok = true;
            return;
        end
    end
end

