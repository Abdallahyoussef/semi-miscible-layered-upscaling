function [stage1_results, layer_ordering, first_param, second_param] = ...
    simulateFirstStage(layer_properties, gamma, ...
                       sw_front, fw_front, slope, SwInitial, mob_total_initial, mob_ratio, ...
                       varargin)
% SIMULATEFIRSTSTAGE Simulate first stage of flooding in stratified layers
%
% DESCRIPTION:
%   This function implements the first stage of unsteady state flooding in stratified
%   porous media according to Youssef & Matthai (2022) methodology.
%   It handles both miscible and immiscible cases and determines layer ordering.
%
% SYNTAX:
%   [stage1_results, layer_ordering] = ...
%       simulateFirstStage(layer_properties, gamma, ...
%                          slope, mob_total_initial, mob_ratio)
%                          
%
%   [stage1_results, layer_ordering] = ...
%       simulateFirstStage(..., 'PropertyName', PropertyValue, ...)
%
% INPUTS:
%   layer_properties    - Structure array with fields:
%       .thickness      - Net pay thickness 
%       .porosity       - Porosity [fraction]
%       .perm_abs       - Absolute permeability
%
%   gamma      - Gamma parameter 
%   slope      - dfw/dsw at front
%   SwInitial - Initial water saturation
%   mob_total_initial   - Total mobility at initial saturation
%   mob_ratio      - Endpoints mobility ratio
%
% OPTIONAL PARAMETERS:
%   'xD_fastest'        - Flood front position in fastest layer (default: 0 to 1 range)
%   'xD_resolution'     - Resolution for xD range (default: 1000)
%   'Verbose'           - Display progress (true/false, default: true)
%   'PlotResults'       - Generate plots (true/false, default: true)
%   'case_name'           - String specifies name for this case
%   'solubility'           - Solubility for miscible cases (mass / volume).
%                                 (default: 0 assuming immiscible flow).
%   'Process'    -  String varaible determing if process is drainage or
%                        imbibition (Drainage, drainage, Imbibition or imbibition, default: imbibition)
%
% OUTPUTS:
%   stage1_results      - Structure with first stage simulation results
%                                   containing list of results 
%                                       average water saturation
%                                       dissolved volume
%                                       free volume
%                                       injected pore volume
%                                       dissolved volume relative to solvent volume
%                                       pore volume injected at end of first stage
%   layer_ordering      - Array with layer indices sorted from fastest to slowest
%
% EXAMPLE:
%   % Get gamma parameters first
%   [gamma_misc, gamma_immisc] = computeGammaParameter(...);
%
%   % Define layer properties
%   layers(1).thickness = 10; layers(1).porosity = 0.25;
%   layers(1).perm_abs = 500; 
%   layers(2).thickness = 8;  layers(2).porosity = 0.22;
%   layers(2).perm_abs = 150; 
%   layers(3).thickness = 3;  layers(3).porosity = 0.15;
%   layers(3).perm_abs = 50;  
%
%   % Simulate first stage
%   [results, ordering] = simulateFirstStage(...
%       layers, gamma_misc, slope_misc, ...
%       mob_total_initial, mob_ratio, 'Verbose', true, 'PlotResults', true);
%
% AUTHOR: AbdAllah AbdRabElNaby Youssef
% DATE: 04 January 2026

% =========================================================================
% PARSE INPUTS
% =========================================================================

% Default parameters
defaults.xD_fastest = [];
defaults.xD_resolution = 1000;
defaults.Verbose = true;
defaults.PlotResults = true;
defaults.case_name = "";
defaults.solubility = 0;
defaults.Process = 'imbibition';


% Parse name-value pairs
p = inputParser;
addRequired(p, 'layer_properties', @isstruct);
addRequired(p, 'gamma', @(x) isnumeric(x) );
addRequired(p, 'slope', @(x) isnumeric(x) );
addRequired(p, 'mob_total_initial', @(x) isnumeric(x) );
addRequired(p, 'mob_ratio', @(x) isnumeric(x) );


fields = fieldnames(defaults);
for i = 1:length(fields)
    addParameter(p, fields{i}, defaults.(fields{i}));
end

parse(p, layer_properties, gamma, slope,  ...
      mob_total_initial, mob_ratio, varargin{:});
      
params = p.Results;

% Validate layer properties
n_layers = length(layer_properties);
validateLayerProperties(layer_properties, n_layers);

% Validate input parameters

% Check if xD_fastest is valid
if any(params.xD_fastest < 0) || any(params.xD_fastest > 1)
    error('simulateFirstStage:UnvalidInput', ...
        'xD_fastest must be postive value between 0 and 1');
end

% Check if resolution is valid
n_points = params.xD_resolution;
if ~isPositiveIntegerArray(n_points)
    error('simulateFirstStage:UnvalidInput', ...
        'xD_resolution must be postive integer');
end

% Check if case name is valid
case_name = params.case_name;
if ~isstring(case_name)
    error('simulateFirstStage:UnvalidInput', ...
        'case name must be a string. Notice that it needs to be between double clotes (e.g., "Case A")');
end

% Check if solubility is valid
solubility = params.solubility;
if ~isnumeric(solubility)
        error('simulateFirstStage:UnvalidInput', ...
        'solubility must be scalar');
elseif solubility < 0     % immiscible flow
    solubility = 0;
end

% Check if process is valid
Process = params.Process;
if ~ischar(Process)
    error('simulateFirstStage:UnvalidInput', ...
        '"Process" must be "Drainage" or "Imbibition"');
elseif strcmp(Process, 'Drainage') || strcmp(Process, 'drainage')
    Process = 'drainage';
elseif strcmp(Process, 'Imbibition') || strcmp(Process, 'imbibition')
    Process = 'imbibition';
else
    error('simulateFirstStage:UnvalidInput', ...
        '"Process" must be "Drainage" or "Imbibition"');
end

% Validate total mobility at initial saturation
if any(mob_total_initial <= 0)
    error('simulateSecondStage:InvalidTotalMobilityAtInitialSaturation', ...
        'Total mobility at initial saturation must be more than 0');
end

if params.Verbose
    fprintf('\n=== First Stage Simulation\n');
    fprintf('gamma: %.6f\n', gamma);
    fprintf('Tangent slope: %.6f\n', slope);
    fprintf('Total mobility at initial saturation: %.6e\n', mob_total_initial);
    fprintf('Number of layers: %d\n', n_layers);
end

% =========================================================================
% COMPUTE FIRST STAGE PARAMETERS
% =========================================================================

if params.Verbose
    fprintf('\n--- Computing First Stage Parameters ---\n');
end

porosity =  [layer_properties.porosity]';
perm_abs = [layer_properties.perm_abs]';
thickness = [layer_properties.thickness]';

% Compute first parameter: γ / (λ_T0 - γ)
first_param = gamma ./ (mob_total_initial - gamma);  % γ / (λ_T0 - γ)
% Compute second parameter: k * f' * λ_T0 / φ
second_param = perm_abs .* mob_total_initial .* slope ./ porosity;      % k * f' * λ_T0 / φ


if params.Verbose
    for i = 1:n_layers
        fprintf('Layer %d:\n', i);
        fprintf('  Thickness: %.2f m, Porosity: %.3f\n', thickness(i), porosity(i));
        fprintf('  Permeability: %.1f mD\n', perm_abs(i));
    end
end


% =========================================================================
% DETERMINE LAYER ORDERING
% =========================================================================

if params.Verbose
    fprintf('\n--- Determining Layer Ordering ---\n');
end

layer_ordering = determineLayerOrdering(first_param, second_param, mob_ratio,   ...
                                       n_layers, params.Verbose);

% Get fastest layer index
fastest_layer_idx = layer_ordering(1);

if params.Verbose
    fprintf('Fastest layer: Layer %d\n', fastest_layer_idx);
    fprintf('Layer ordering (fastest to slowest):\n');
    for i = 1:n_layers
        fprintf('  %d: Layer %d\n', i, layer_ordering(i));
    end
end

% =========================================================================
% COMPUTE FLOOD FRONT ADVANCEMENT
% =========================================================================

if params.Verbose
    fprintf('\n--- Computing Flood Front Advancement ---\n');
end

% If xD_fastest is provided, use it; otherwise use full range
if isempty(params.xD_fastest)
    xD_fastest_range = linspace(0, 1, n_points)';
else
    xD_fastest_range = params.xD_fastest(:);
end

% Compute flood front positions in all layers
xD_all_layers = computeFrontPosition(xD_fastest_range, first_param, second_param, ...
                                           first_param(fastest_layer_idx), second_param(fastest_layer_idx), mob_ratio);


% =========================================================================
% COMPUTE AVERAGE SATURATIONS AND VOLUMES
% =========================================================================

if params.Verbose
    fprintf('\n--- Computing Average Saturations and Volumes of First Stage ---\n');
end

[avg_saturation, free_injected_volume, soluble_injected_volume, injected_pv, solute_solvent_volume_ratio] = ...
    computeVolumes( xD_all_layers, sw_front, fw_front, slope, ...
   porosity, thickness,  SwInitial, solubility, ...
    Process, params.Verbose);

% =========================================================================
% COMPUTE FLOW RATE PARTITIONING & PRESSURE DROP
% =========================================================================
% Compute effective total mobility
mob_total_effective_reciprocal = (xD_all_layers ./ gamma) + (1 ./ mob_total_initial) .* (1 - xD_all_layers);
mob_total_effective = 1 ./ mob_total_effective_reciprocal;
% Compute layer flow capacities
flow_capacity = (perm_abs .* thickness) .* mob_total_effective;
% Compute flow rate ratio relative to fastest layer
layer_rate_realtive_to_fastest_layer = flow_capacity ./ flow_capacity(fastest_layer_idx, :);
% Compute flow rate of fastest layer
fasest_layer_rate = 1 ./sum(layer_rate_realtive_to_fastest_layer);
% Compute flow rate in all layers
rate_layers = fasest_layer_rate .* layer_rate_realtive_to_fastest_layer;
% Compute pressure drop
pressure_drop = fasest_layer_rate ./ flow_capacity(fastest_layer_idx, :);

% =========================================================================
% ORGANIZE RESULTS
% =========================================================================

stage1_results = struct();
stage1_results.case_name = case_name;
stage1_results.solubility = solubility;
stage1_results.gamma_used = gamma;
stage1_results.slope_used = slope;
stage1_results.mob_total_initial = mob_total_initial;

% Store parameters
stage1_results.first_param = first_param;
stage1_results.second_param = second_param;

% Store layer information
stage1_results.layer_properties = layer_properties;
stage1_results.layer_ordering = layer_ordering;
stage1_results.fastest_layer_idx = fastest_layer_idx;

% Store flood front advancement
stage1_results.xD_fastest_range = xD_fastest_range;
stage1_results.xD_all_layers = xD_all_layers;

% Store computed masses and saturations
stage1_results.avg_saturation = avg_saturation;
stage1_results.free_injected_volume = free_injected_volume;
stage1_results.soluble_injected_volume = soluble_injected_volume;
stage1_results.injected_pv = injected_pv;
stage1_results.solute_solvent_volume_ratio = solute_solvent_volume_ratio;

% Store computed pressure drop and flow rate and saturations
stage1_results.pressure_drop = pressure_drop;
stage1_results.layer_rate = rate_layers;

% =========================================================================
% GENERATE PLOTS
% =========================================================================

if params.PlotResults
    visualizeFirstStageResults(stage1_results, params.Verbose);
end

if params.Verbose
    fprintf('\n=== First Stage Simulation Complete ===\n\n');
    displayStage1Summary(stage1_results);
end

end % Main function

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function validateLayerProperties(layers, n_layers)
% VALIDATELAYERPROPERTIES Validate layer property structure
    
    required_fields = {'thickness', 'porosity', 'perm_abs'};
    
    for i = 1:n_layers
        for j = 1:length(required_fields)
            if ~isfield(layers(i), required_fields{j})
                error('simulateFirstStage:MissingField', ...
                      'Layer %d missing field: %s', i, required_fields{j});
            end
        end
        
        % Validate values
        if layers(i).thickness <= 0
            error('simulateFirstStage:InvalidThickness', ...
                  'Layer %d: thickness must be positive', i);
        end
        
        if layers(i).porosity <= 0 || layers(i).porosity >= 1
            error('simulateFirstStage:InvalidPorosity', ...
                  'Layer %d: porosity must be between 0 and 1', i);
        end
        
        if layers(i).perm_abs <= 0
            error('simulateFirstStage:InvalidPermeability', ...
                  'Layer %d: absolute permeability must be positive', i);
        end
    end
end

function xD_layer_i = computeFrontPosition(xD_fastest, first_param_i, second_param_i, ...
                                           first_param_fastest, second_param_fastest, mob_ratio_i)
% COMPUTEFRONTPOSITION Compute front position in layer i using Eq. from paper
%   Vectorized version that computes for all combinations of xD_fastest values and layers
    
    % Extract parameters
    first_param_i = first_param_i(:);  % Ensure column vector [n_layers × 1]
    second_param_i = second_param_i(:);  
    xD_fastest = xD_fastest(:)';  % Ensure row vector [1 × n_xD_values]
    

    % Compute the term inside square root (using implicit expansion)
    term1 = (0.5 / first_param_fastest) *  (xD_fastest .^2) + xD_fastest;  
    term2 =(2 ./  first_param_i).* second_param_i / second_param_fastest;

   
    radicand = 1 + term1 .* term2;

     
    % Check for negative radicands
    negative_mask = radicand < 0;
    if any(negative_mask(:))
        warning('simulateFirstStage:NegativeRadicand', ...
                'Negative radicand detected in %d combinations. Using absolute value.', ...
                sum(negative_mask(:)));
        radicand(negative_mask) = abs(radicand(negative_mask));
    end

    % Choosing correct formula of position in layer i based on endpoint mobility ratio
    if all(mob_ratio_i ~= 1)
        xD_layer_i = first_param_i .* (-1 + sqrt(radicand));     
    elseif all(mob_ratio_i == 1)
        term2 = second_param_i / second_param_fastest;
        xD_layer_i = term1 .* term2;
    else
        n_layers = length(first_param_i);
        n_xD_fastest = length(xD_fastest);
        xD_layer_i = zeros(n_layers , n_xD_fastest);
        for i = 1:n_layers
            if mob_ratio_i(i) ~= 1
                xD_layer_i(i,:) = first_param_i(i) .* (-1 + sqrt(radicand(i,:)));   
            else
                term2(i) = second_param_i(i) / second_param_fastest;
                xD_layer_i(i,:) = term1 .* term2(i);
            end
        end
    end
    
    % Checking if there are some vlaues of flood front position are
    % non-physical (negative)
    % if any(xD_layer_i < 0)
    %     [non_physical_layer_idx, ~] = find(xD_layer_i < 0);
    %     non_physical_layer_idx = unique(non_physical_layer_idx);
    %     error('simulateFirstStage:NonPhysicalValues', ...
    %                   'Layer %d gives non physical flood front position', non_physical_layer_idx);
    % end

end


function ordering = determineLayerOrdering(first_param, second_param, ...
                                          mob_ratio, n_layers, verbose)
% DETERMINELAYERORDERING Determine layer ordering from fastest to slowest
    
    if verbose
        fprintf('  Testing all possible layer orderings...\n');
    end

    
    % Initialize maximum front position
    max_xD_test = inf;
    fastest_idx = 0;


    while max_xD_test > 1  && fastest_idx < n_layers
        fastest_idx = fastest_idx + 1;
        first_param_fastest = first_param(fastest_idx);
        second_param_fastest = second_param(fastest_idx);

        % Test with xD_fastest = 1
        xD_fastest = 1;
        xD_test = computeFrontPosition(xD_fastest, first_param, second_param, ...
                                           first_param_fastest, second_param_fastest, mob_ratio);
        xD_test(fastest_idx) = 1;
        
        % Sorting layers from fastest to slowest
         [~,ordering] = sort(xD_test, 'descend');
        
         % Find maximum xD (should be 1 for fastest, others < 1)
          max_xD_test = max(xD_test);
    end
end

function [avg_sat_total, free_phase, soluble_phase, injected_pv, solute_solvent_volume_ratio] = computeVolumes(...
    xD_all_layers, sw_front, fw_front, slope, porosity, thickness,  ...
    sw_initial, solubility, Process, verbose)
% COMPUTEVolumes Compute average saturations and volumes dissoloved and
% free

    
% Compute pore volumes 
    pore_volume = porosity .* thickness;     % normalized pore volume of each layer
    total_pore_volume = sum(pore_volume);  % total pore volume
    
    if verbose
        fprintf('Total pore volume (per unit width and unit length): %.4f L³/L²\n', total_pore_volume);
    end
    
% Compute injected volumes
    injected_volume = (thickness .* porosity ./ slope) .* xD_all_layers;     % This estimation is correct for immiscible and approximation for immiscible 
    total_injecetd_volume = sum(injected_volume);       

% Compute pore volume injected
    injected_pv = total_injecetd_volume / total_pore_volume;


% Compute average saturation for each layer
    if strcmp(Process, 'drainage')
        fw_inlet = 0;
    else
        fw_inlet = 1;
    end
    sw_avg_behind = sw_front - (1 ./ slope) .* (fw_front - fw_inlet);     % average saturation behind front
    sw_avg = sw_avg_behind .* xD_all_layers + sw_initial .* (1 - xD_all_layers);   % average saturation in each layer

% Compute water volume
    water_volume = (thickness .* porosity) .* sw_avg;
    water_volume_total = sum(water_volume);

% Compute model avearge saturation of whole model
    avg_sat_total = water_volume_total / total_pore_volume;

% Determine free phase volume
   if strcmp(Process, 'drainage')
       free_phase = 1 - avg_sat_total;
   else
       free_phase = avg_sat_total;
   end

% Determine if case miscble or no
    if solubility > 0
        soluble_phase = injected_pv - free_phase;
        if strcmp(Process, 'drainage')
            solute_solvent_volume_ratio = ((1 ./ slope) - (1 - sw_avg_behind)) ./ sw_avg_behind;
        else
            solute_solvent_volume_ratio = ((1 ./ slope) - sw_avg_behind) ./ (1 - sw_avg_behind);
        end
        
    else 
        soluble_phase = 0 * injected_pv;
        solute_solvent_volume_ratio = 0 * porosity;
    end
    


    if verbose
        fprintf('  Final pore volume of injected component free at end of first stage: %.4f \n', free_phase(end));
        fprintf('  Final pore volume of injected component soluble at end of first stage: %.4f \n', soluble_phase(end));
        fprintf('  Solute to solvent volume rattio: %.4f\n', solute_solvent_volume_ratio(end));
        fprintf('  Final injected pore volume at end of first stage: %.4f\n', injected_pv(end));
    end
end


function visualizeFirstStageResults(results, verbose)
% VISUALIZEFIRSTSTAGERESULTS Create comprehensive plots for first stage
    
    if verbose
        fprintf('Generating first stage plots...\n');
    end
    
    fig = figure('Position', [50, 50, 1400, 900], 'Color', 'white');
    
    % Colors for different layers
    colors = lines(results.layer_ordering(end));
    
    % Subplot 1: Flood front advancement
    subplot(2, 3, 1);
    hold on;
    
    for i = 1:length(results.layer_ordering)
        layer_idx = results.layer_ordering(i);
        plot(results.xD_fastest_range, results.xD_all_layers(layer_idx, :), ...
             '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Layer %d', layer_idx));
    end
    
    plot([0, 1], [0, 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Reference');
    
    xlabel('xD_{fastest}');
    ylabel('xD_i');
    title('Flood Front Advancement');
    legend('Location', 'best', 'NumColumns', 2);
    grid on;
    box on;
    xlim([0, 1]);
    ylim([0, 1]);
    
    % Subplot 2: Free volume vs injected PV
    subplot(2, 3, 2);
    plot(results.injected_pv, results.free_injected_volume, 'b-', 'LineWidth', 2);
    xlabel('Injected Pore Volume');
    ylabel('Free Volume ');
    title('Free Volume vs Injection');
    grid on;
    box on;
    
    % Subplot 3: Average saturation profile
    subplot(2, 3, 3);
    
    plot(results.xD_fastest_range, results.avg_saturation, 'r-', 'LineWidth', 2);
    
    xlabel('xD_{fastest}');
    ylabel('<s_w>');
    title('Domain Average Saturation');
    grid on;
    box on;
    xlim([0, 1]);

    
    % Subplot 4: Pressure drop
    subplot(2, 3, 4);
    
    plot(results.injected_pv, results.pressure_drop, 'b-', 'LineWidth', 2);
    xlabel('Injected Pore Volume');
    ylabel('\deltaP');
    title('Pressure Drop');
    grid on;
    box on;
    

    % Subplot 5: Flood front advancement vs injected pore volume
    subplot(2, 3, 5);
    hold on;
    
    for i = 1:length(results.layer_ordering)
        layer_idx = results.layer_ordering(i);
        plot(results.injected_pv, results.xD_all_layers(layer_idx, :), ...
             '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Layer %d', layer_idx));
    end
    
    xlabel('Injected Pore Volume');
    ylabel('xD_i');
    title('Flood Front Advancement');
    legend('Location', 'best', 'NumColumns', 2);
    grid on;
    box on;
    ylim([0, 1]);

    % Subplot 6: dissolved volume vs injected PV
    subplot(2, 3, 6);
    plot(results.injected_pv, results.soluble_injected_volume, 'g-', 'LineWidth', 2);
    xlabel('Injected Pore Volume');
    ylabel('Dissolved Volume ');
    title('Dissolved Volume vs Injection');
    grid on;
    box on;

   
    sgtitle(sprintf('First Stage: %s', results.case_name), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    if verbose
        fprintf('  Plots generated successfully.\n');
    end
end

function displayStage1Summary(results)
% DISPLAYSTAGE1SUMMARY Display summary of first stage results
    
    fprintf('\n=== First Stage Simulation Summary ===\n');
    fprintf('Case: %s\n', results.case_name);
    fprintf('Gamma parameter used: %.6f\n', results.gamma_used);
    fprintf('Fractional flow slope used: %.4f\n', results.slope_used);
    fprintf('Total mobility at initial saturation: %.6e\n', results.mob_total_initial);
    fprintf('\nLayer ordering (fastest to slowest):\n');
    
    for i = 1:length(results.layer_ordering)
        fprintf('  %d. Layer %d\n', i, results.layer_ordering(i));
    end
    
    fprintf('\nFinal results (xD_fastest = 1):\n');
    fprintf('  Average water saturation: %.4f\n', results.avg_saturation(end));
    fprintf('  Free Volume: %.4f \n', results.free_injected_volume(end));
    fprintf('  Soluble Volume: %.2f \n', results.soluble_injected_volume(end));
    
end