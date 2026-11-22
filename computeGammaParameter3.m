function [gamma_miscible, gamma_immiscible, integration_details] = ...
    computeGammaParameter3(Swr, Snwr, sw_front_miscible, sw_front_immiscible, ...
                          fw_func, dfw_dsw_func, mob_t_func, varargin)
% COMPUTEGAMMAPARAMETER Compute gamma parameter for miscible and immiscible cases
%
% DESCRIPTION:
%   This function computes the gamma parameter as defined in the Youssef & Matthai
%   (2022) paper for CO2 flooding. Gamma is defined as the integral from Swr to
%   the flood front saturation of the second derivative of fractional flow
%   divided by total mobility. The function now supports both single and multiple
%   layers with comprehensive visualization.
%
% IMPORTANT NOTE:
%   For drainage process: Integration is from Swr to sw_front
%   For imbibition process: Integration is from 1-Snwr to sw_front
%
% SYNTAX:
%   For single layer:
%   [gamma_miscible, gamma_immiscible] = ...
%       computeGammaParameter(Swr, Snwr, sw_front_miscible, sw_front_immiscible, ...
%                             fw_func, dfw_dsw_func, mob_total_func)
%
%   For multiple layers:
%   [gamma_miscible, gamma_immiscible] = ...
%       computeGammaParameter(Swr, Snwr, sw_front_miscible, sw_front_immiscible, ...
%                             fw_func, dfw_dsw_func, mob_total_func, layer_idx)
%
%   [gamma_miscible, gamma_immiscible, integration_details] = ...
%       computeGammaParameter(..., 'PropertyName', PropertyValue, ...)
%
% INPUTS:
%   Swr                 - Residual water saturation (scalar or vector for multiple layers)
%   Snwr                - Residual non-wetting phase saturation (scalar or vector)
%   sw_front_miscible   - Flood front saturation for miscible case (scalar or vector)
%   sw_front_immiscible - Flood front saturation for immiscible case (scalar or vector)
%   fw_func             - Function handle or cell array: fw(sw)
%   dfw_dsw_func        - Function handle or cell array: dfw_dsw(sw)
%   mob_t_func          - Function handle or cell array: λ_total(sw)
%
% OPTIONAL INPUTS:
%   layer_idx    - (Optional) Layer index for multi-layer systems (default: all layers)
%
% OPTIONAL PARAMETERS:
%   'IntegrationPoints' - Number of points for numerical integration (default: 1000)
%   'IntegrationMethod' - 'trapz' or 'quad' (default: 'quad')
%   'SecondDerivativeMethod' -  'complex', or 'finite' (default: 'finite')
%   'ComplexStepSize'   - Step size for complex step derivative (default: 1e-20)
%   'FiniteDiffStep'    - Step size for finite difference (default: 1e-6)
%   'Verbose'           - Display progress and details (default: true)
%   'PlotIntegrand'     - Plot the integrand function (default: true)
%   'CombinePlots'      - Combine multiple layers in single figure with subplots (default: true)
%   'LayerNames'        - Cell array of layer names for labeling (default: {})
%   'Compute'           - Which gamma values to compute: 'both', 'miscible', or 'immiscible' (default: 'both')
%   'Process'           - String varaible determing if process is drainage or
%                         imbibition (Drainage, drainage, Imbibition or imbibition, default: imbibition)
%
% OUTPUTS:
%   gamma_miscible      - Gamma parameter for miscible case (scalar or vector)
%   gamma_immiscible    - Gamma parameter for immiscible case (scalar or vector)
%   integration_details - Structure with detailed integration results
%
% EXAMPLE:
%   % Drainage process (miscible only)
%   [gamma_misc, ~] = computeGammaParameter(...
%       params.Swr, params.Snwr, sw_front_miscible, [], ...
%       fw_func, dfw_dsw_func, mob_t_func, ...
%       'Compute', 'miscible', 'Process', 'drainage');
%
%   % Imbibition process (both)
%   [gamma_misc, gamma_immisc] = computeGammaParameter(...
%       params.Swr, params.Snwr, sw_front_miscible, sw_front_immiscible, ...
%       fw_func, dfw_dsw_func, mob_t_func, ...
%       'Process', 'imbibition');
%
% AUTHOR: AbdAllah AbdRabElnaby Youssef ElSaieed
% DATE: 2026-01-05
% MODIFIED: To support multiple layers, selective computation, and process type (2026-01-14)

% =========================================================================
% PARSE INPUTS
% =========================================================================

% Default parameters
defaults.IntegrationPoints = 1000;
defaults.IntegrationMethod = 'quad';
defaults.SecondDerivativeMethod = 'finite';
defaults.ComplexStepSize = 1e-20;
defaults.FiniteDiffStep = 1e-6;
defaults.Verbose = true;
defaults.PlotIntegrand = true;
defaults.CombinePlots = true;
defaults.LayerNames = {};
defaults.Compute = 'both';  % 'both', 'miscible', or 'immiscible'
defaults.Process = 'imbibition';  % 'drainage' or 'imbibition'

% Check if layer_idx is provided as numeric argument
layer_idx = [];      % initializing layer id with empty vector
if ~isempty(varargin) &&  isscalar(varargin{1}) &&  isPositiveIntegerArray(varargin{1})
    layer_idx = varargin{1};
    varargin = varargin(2:end);
end

% Parse name-value pairs
p = inputParser;
addRequired(p, 'Swr', @(x) isnumeric(x));
addRequired(p, 'Snwr', @(x) isnumeric(x));
addRequired(p, 'sw_front_miscible', @(x) isnumeric(x) || isempty(x));
addRequired(p, 'sw_front_immiscible', @(x) isnumeric(x) || isempty(x));
addRequired(p, 'fw_func');
addRequired(p, 'dfw_dsw_func');
addRequired(p, 'mob_t_func');

fields = fieldnames(defaults);
for i = 1:length(fields)
    addParameter(p, fields{i}, defaults.(fields{i}));
end

parse(p, Swr, Snwr, sw_front_miscible, sw_front_immiscible, ...
      fw_func, dfw_dsw_func, mob_t_func, varargin{:});
params = p.Results;

% Validate Compute parameter
valid_compute = {'both', 'miscible', 'immiscible'};
params.Compute = validatestring(params.Compute, valid_compute);

% Validate Process parameter
valid_process = {'drainage', 'imbibition', 'Drainage' , 'Imbibition'};
params.Process = validatestring(params.Process, valid_process);

if strcmp(params.Process, 'Drainage') || strcmp(params.Process, 'drainage')
    warning('findShockFront:Process', ...
        'This is draiange process. If you need to switch to imbibition, set "Process" to "imbibition"');

    params.Process = 'drainage';

elseif strcmp(params.Process, 'Imbibition') || strcmp(params.Process, 'imbibition')
    warning('findShockFront:Process', ...
        'This is imbibition process. If you need to switch to draiange, set "Process" to "drainage"');

        params.Process = 'imbibition';
end
   


% Validate that appropriate inputs are provided based on Compute parameter
if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
    if isempty(sw_front_miscible)
        error('computeGammaParameter:MissingInput', ...
              'sw_front_miscible must be provided when Compute is ''%s''', params.Compute);
    end
end

if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
    if isempty(sw_front_immiscible)
        error('computeGammaParameter:MissingInput', ...
              'sw_front_immiscible must be provided when Compute is ''%s''', params.Compute);
    end
end

% Validate integration method
valid_methods = {'trapz', 'quad'};
params.IntegrationMethod = validatestring(params.IntegrationMethod, valid_methods);

% Validate second derivative method
valid_deriv_methods = {'complex', 'finite'};
params.SecondDerivativeMethod = validatestring(params.SecondDerivativeMethod, ...
                                              valid_deriv_methods);

% =========================================================================
% DETERMINE NUMBER OF LAYERS, SETUP AND QC INPUTS
% =========================================================================

% Check if we have multiple layers
if iscell(fw_func)
    n_layers = length(fw_func);
    
    % Check if the calculation of all layers are required
    if isempty(layer_idx)
        % Validate that derivative functions match
        if ~iscell(dfw_dsw_func) || length(dfw_dsw_func) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of fractional flow functions must match derivative functions');

        % Validate that mobility functions match
        elseif ~iscell(mob_t_func) || length(mob_t_func) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of mobility functions must match fractional flow functions');

        % Validate that residual water saturations match
        elseif ~isvector(Swr) || length(Swr) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of residual water saturations must match number of layers');

        % Validate that residual non-wetting saturations match
        elseif ~isvector(Snwr) || length(Snwr) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of residual non-wetting saturations must match number of layers');

        end
    
    % Check if calculations are required for specific layers
    elseif ~isPositiveIntegerArray(layer_idx) || max(layer_idx) > n_layers
        error('computeGammaParameter:InvalidInput', ...
            'layer_idx must be positive integer and less than number of layers')

    else   % only specific layers are required
        n_layers = length(layer_idx);

        % Validate that residual saturations match
        if ~isvector(Swr) || length(Swr) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of residual water saturations must match layer_idx');
        end
        
        if ~isvector(Snwr) || length(Snwr) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of residual non-wetting saturations must match layer_idx');
        end

        % Extract specific layer functions
        fw_func = fw_func(layer_idx);
        dfw_dsw_func = dfw_dsw_func(layer_idx);
        mob_t_func = mob_t_func(layer_idx);
    end

    % Validate miscible front saturations if needed
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
        if ~isvector(sw_front_miscible) || length(sw_front_miscible) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of miscible front saturations must match number of layers');
        end
    end
    
    % Validate immiscible front saturations if needed
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
        if ~isvector(sw_front_immiscible) || length(sw_front_immiscible) ~= n_layers
            error('computeGammaParameter:InvalidInput', ...
                  'Number of immiscible front saturations must match number of layers');
        end
    end

else
    n_layers = 1;
    
    % Convert to cell arrays for uniform handling
    fw_func = {fw_func};
    dfw_dsw_func = {dfw_dsw_func};
    mob_t_func = {mob_t_func};
    
    % Ensure saturations are vectors
    if isscalar(Swr)
        Swr = Swr(:);
    end
    if isscalar(Snwr)
        Snwr = Snwr(:);
    end
    if ~isempty(sw_front_miscible) && isscalar(sw_front_miscible)
        sw_front_miscible = sw_front_miscible(:);
    end
    if ~isempty(sw_front_immiscible) && isscalar(sw_front_immiscible)
        sw_front_immiscible = sw_front_immiscible(:);
    end
end

% Calculate lower bound based on process
if strcmp(params.Process, 'drainage')
    lower_bound = Swr;  % For drainage: integrate from Swr
    lower_bound_label = 'Swr';
else
    lower_bound = 1 - Snwr;  % For imbibition: integrate from 1-Snwr
    lower_bound_label = '1-Snwr';
end

% Validate saturation bounds for all layers based on what's being computed
for layer = 1:n_layers
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
        % Validate against appropriate lower bound
        if sw_front_miscible(layer) >= 1 - Snwr(layer) || sw_front_miscible(layer) <= Swr(layer)
            error('computeGammaParameter:InvalidFrontMiscible', ...
                  'Layer %d: Miscible front saturation must be in (%s, %s), got %.4f', ...
                  layer, Swr(layer), 1 - Snwr(layer), sw_front_miscible(layer));
        end
    end
    
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
        if sw_front_immiscible(layer) >= 1 - Snwr(layer) || sw_front_immiscible(layer) <= Swr(layer)
            error('computeGammaParameter:InvalidFrontImmiscible', ...
                  'Layer %d: Immiscible front saturation must be in (%s, %s), got %.4f', ...
                  layer, Swr(layer), 1 - Snwr(layer), sw_front_immiscible(layer));
        end
    end
end

% Check if input endpoint saturations are valid
if min(Swr) < 0 || max(Swr) >= 1 || min(Snwr) < 0 || max(Snwr) >= 1 || max(Swr + Snwr) >= 1
        error('computeDoubleIntegralTerm:InvalidResiduals', ...
              ['Invalid residual saturations. Input residual saturations must be between zero and 1, ' ...
              '& their sum less than or equal 1']);
end

if params.Verbose
    fprintf('\n=== Computing Gamma Parameter ===\n');
    fprintf('Number of layers: %d\n', n_layers);
    fprintf('Process: %s\n', params.Process);
    fprintf('Integration from: %s\n', lower_bound_label);
    fprintf('Computation mode: %s\n', params.Compute);
    fprintf('Integration method: %s\n', params.IntegrationMethod);
    fprintf('Second derivative method: %s\n', params.SecondDerivativeMethod);
    if n_layers > 1
        fprintf('Lower bounds (%s): %s\n', lower_bound_label, mat2str(lower_bound', 3));
    else
        fprintf('Lower bound (%s): %.4f\n', lower_bound_label, lower_bound(1));
    end
end


% =========================================================================
% PROCESS EACH LAYER
% =========================================================================

% Initialize outputs based on what's being computed
if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
    gamma_miscible = zeros(n_layers, 1);
else
    gamma_miscible = [];
end

if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
    gamma_immiscible = zeros(n_layers, 1);
else
    gamma_immiscible = [];
end

% Prepare results structure
integration_details = struct();
integration_details.n_layers = n_layers;
integration_details.Compute = params.Compute;
integration_details.Process = params.Process;
integration_details.Swr = Swr;
integration_details.Snwr = Snwr;
integration_details.lower_bound = lower_bound;
integration_details.lower_bound_label = lower_bound_label;

% Initialize based on what's being computed
if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
    integration_details.gamma_miscible = zeros(n_layers, 1);
    integration_details.sw_front_miscible = sw_front_miscible;
    integration_details.details_miscible = cell(n_layers, 1);
end

if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
    integration_details.gamma_immiscible = zeros(n_layers, 1);
    integration_details.sw_front_immiscible = sw_front_immiscible;
    integration_details.details_immiscible = cell(n_layers, 1);
end

integration_details.methods = struct(...
    'integration_method', params.IntegrationMethod, ...
    'second_derivative_method', params.SecondDerivativeMethod);

% Looping on layers to compute gamma parameters
for layer = 1:n_layers
    if params.Verbose
        fprintf('\n--- Processing Layer %d ---\n', layer);
        fprintf('Process: %s\n', params.Process);
        fprintf('Integration from %s: %.4f', lower_bound_label, lower_bound(layer));
        if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
            fprintf(', Miscible front: %.4f', sw_front_miscible(layer));
        end
        if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
            fprintf(', Immiscible front: %.4f', sw_front_immiscible(layer));
        end
        fprintf('\n');
    end
    
    % Get functions for this layer
    current_fw_func = @(sw) fw_func{layer}(sw);
    current_dfw_func = @(sw) dfw_dsw_func{layer}(sw);
    current_mob_t_func = @(sw) mob_t_func{layer}(sw);
    
    % =====================================================================
    % CREATE SECOND DERIVATIVE FUNCTION
    % =====================================================================
    
    if params.Verbose
        fprintf('Creating second derivative function...\n');
    end
    
    % Create function for second derivative based on selected method
    d2fw_dsw2_func = createSecondDerivativeFunction(...
        current_fw_func, params.SecondDerivativeMethod, ...
        params.ComplexStepSize, params.FiniteDiffStep, params.Verbose);
    
    % =====================================================================
    % CREATE INTEGRAND FUNCTION
    % =====================================================================
    
    if params.Verbose
        fprintf('Creating integrand function: f''''(s) / λ_T(s)\n');
    end
    
    % Integrand function: f''(s) / λ_T(s)
    integrand_func = @(sw) computeIntegrand(sw, d2fw_dsw2_func, current_mob_t_func);
    
    % =====================================================================
    % COMPUTE GAMMA FOR MISCIBLE CASE (if needed)
    % =====================================================================
    
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible')
        if params.Verbose
            fprintf('Computing gamma for miscible case...\n');
            fprintf('Integration from %.4f to %.4f\n', lower_bound(layer), sw_front_miscible(layer));
        end
        
        [gamma_int_misc, details_miscible] = computeGammaIntegral(...
            lower_bound(layer), sw_front_miscible(layer), integrand_func, params.IntegrationPoints, ...
             params.IntegrationMethod, params.Verbose);
        
        gamma_miscible_layer = current_dfw_func(sw_front_miscible(layer)) ./ gamma_int_misc;
        
        % Store results
        gamma_miscible(layer) = gamma_miscible_layer;
        integration_details.gamma_miscible(layer) = gamma_miscible_layer;
        integration_details.details_miscible{layer} = details_miscible;
        integration_details.integrand_func{layer} = integrand_func;
        
        if params.Verbose
            fprintf('  Miscible gamma: %.6f\n', gamma_miscible_layer);
        end
    end
    
    % =====================================================================
    % COMPUTE GAMMA FOR IMMISCIBLE CASE (if needed)
    % =====================================================================
    
    if strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible')
        if params.Verbose
            fprintf('Computing gamma for immiscible case...\n');
            fprintf('Integration from %.4f to %.4f\n', lower_bound(layer), sw_front_immiscible(layer));
        end
        
        [gamma_int_immisc, details_immiscible] = computeGammaIntegral(...
            lower_bound(layer), sw_front_immiscible(layer), integrand_func, params.IntegrationPoints, ...
             params.IntegrationMethod, params.Verbose);
        
        gamma_immiscible_layer = current_dfw_func(sw_front_immiscible(layer)) ./ gamma_int_immisc;
        
        % Store results
        gamma_immiscible(layer) = gamma_immiscible_layer;
        integration_details.gamma_immiscible(layer) = gamma_immiscible_layer;
        integration_details.details_immiscible{layer} = details_immiscible;
        if ~exist('integrand_func', 'var')  % Only store if not already stored
            integration_details.integrand_func{layer} = integrand_func;
        end
        
        if params.Verbose
            fprintf('  Immiscible gamma: %.6f\n', gamma_immiscible_layer);
        end
    end
    
    % Display ratio only if both are computed
    if strcmp(params.Compute, 'both')
        if params.Verbose
            ratio = gamma_miscible(layer) / gamma_immiscible(layer);
            fprintf('  Ratio (M/I): %.4f\n', ratio);
        end
    end
end

% =========================================================================
% DISPLAY FINAL RESULTS
% =========================================================================

if params.Verbose
    fprintf('\n=== Gamma Parameter Results ===\n');
    fprintf('Process: %s\n', params.Process);
    
    switch params.Compute
        case 'both'
            if n_layers == 1
                fprintf('Miscible case: γ = %.6f\n', gamma_miscible(1));
                fprintf('Immiscible case: γ = %.6f\n', gamma_immiscible(1));
                fprintf('Ratio (Miscible/Immiscible): %.4f\n', ...
                        gamma_miscible(1) / gamma_immiscible(1));
                fprintf('Difference: %.6f\n', gamma_miscible(1) - gamma_immiscible(1));
            else
                fprintf('Layer\tγ_miscible\tγ_immiscible\tRatio (M/I)\n');
                fprintf('-----\t----------\t------------\t-----------\n');
                for layer = 1:n_layers
                    ratio = gamma_miscible(layer) / gamma_immiscible(layer);
                    fprintf('%d\t%.6f\t%.6f\t%.4f\n', ...
                            layer, gamma_miscible(layer), gamma_immiscible(layer), ratio);
                end
                
                fprintf('\nSummary Statistics:\n');
                fprintf('  Average γ_miscible: %.6f\n', mean(gamma_miscible));
                fprintf('  Average γ_immiscible: %.6f\n', mean(gamma_immiscible));
                fprintf('  Average ratio (M/I): %.4f\n', mean(gamma_miscible ./ gamma_immiscible));
            end
            
        case 'miscible'
            if n_layers == 1
                fprintf('Miscible gamma: γ = %.6f\n', gamma_miscible(1));
            else
                fprintf('Layer\tγ_miscible\n');
                fprintf('-----\t----------\n');
                for layer = 1:n_layers
                    fprintf('%d\t%.6f\n', layer, gamma_miscible(layer));
                end
                
                fprintf('\nSummary Statistics:\n');
                fprintf('  Average: %.6f\n', mean(gamma_miscible));
                fprintf('  Min: %.6f\n', min(gamma_miscible));
                fprintf('  Max: %.6f\n', max(gamma_miscible));
                fprintf('  Std: %.6f\n', std(gamma_miscible));
            end
            
        case 'immiscible'
            if n_layers == 1
                fprintf('Immiscible gamma: γ = %.6f\n', gamma_immiscible(1));
            else
                fprintf('Layer\tγ_immiscible\n');
                fprintf('-----\t------------\n');
                for layer = 1:n_layers
                    fprintf('%d\t%.6f\n', layer, gamma_immiscible(layer));
                end
                
                fprintf('\nSummary Statistics:\n');
                fprintf('  Average: %.6f\n', mean(gamma_immiscible));
                fprintf('  Min: %.6f\n', min(gamma_immiscible));
                fprintf('  Max: %.6f\n', max(gamma_immiscible));
                fprintf('  Std: %.6f\n', std(gamma_immiscible));
            end
    end
end

% =========================================================================
% PLOT INTEGRAND (only if PlotIntegrand is true and we have data to plot)
% =========================================================================

if params.PlotIntegrand
    % Determine if we have data to plot
    has_miscible_data = strcmp(params.Compute, 'both') || strcmp(params.Compute, 'miscible');
    has_immiscible_data = strcmp(params.Compute, 'both') || strcmp(params.Compute, 'immiscible');
    
    if has_miscible_data || has_immiscible_data
        if params.Verbose
            fprintf('\nGenerating integrand plots...\n');
        end
        plotGammaIntegrandMultiLayer(integration_details, params.Verbose, params.CombinePlots, params.LayerNames);
    else
        if params.Verbose
            fprintf('\nNo data to plot (both miscible and immiscible computations are disabled)\n');
        end
    end
end

if params.Verbose
    fprintf('\n=== Gamma Computation Complete ===\n\n');
end

end % Main function

% =========================================================================
% HELPER FUNCTIONS (unchanged)
% =========================================================================

function d2fw_dsw2_func = createSecondDerivativeFunction(...
    fw_func, method, complex_step, finite_step, verbose)
% CREATESECONDDERIVATIVEFUNCTION Create function for second derivative
    
    switch method
        case 'complex'
            % Complex step method (most accurate)
            if verbose
                fprintf('  Using complex step method (step = %.1e)\n', complex_step);
            end
            
            d2fw_dsw2_func = @(sw) complexStepSecondDerivative(...
                fw_func, sw, complex_step);
            
        case 'finite'
            % Finite difference method
            if verbose
                fprintf('  Using finite difference method (step = %.1e)\n', finite_step);
            end
            
            d2fw_dsw2_func = @(sw) finiteDifferenceSecondDerivative(...
                fw_func, sw, finite_step);
            
        otherwise
            error('computeGammaParameter:InvalidMethod', ...
                  'Unknown second derivative method: %s', method);
    end
end

function d2fw_dsw2 = complexStepSecondDerivative(fw_func, sw, h)
% COMPLEXSTEPSECONDDERIVATIVE Compute second derivative using complex step
    
    % Direct complex step for second derivative
    % f''(x) ≈ ( f(x) - f(x+ih) ) / 2h^2
    f_fun = fw_func(sw);
    f_plus = real(fw_func(sw + h * 1i));
    d2fw_dsw2 = 2 * (f_fun - f_plus) / (h^2);
end

function d2fw_dsw2 = finiteDifferenceSecondDerivative(fw_func, sw, h)
% FINITEDIFFERENCESECONDDERIVATIVE Compute second derivative using finite difference
    
    % Central difference for second derivative
    % f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h^2
    df_plus = fw_func(sw + h);
    df_minus = fw_func(sw - h);
    df_center = fw_func(sw);
    
    d2fw_dsw2 = (df_plus - 2*df_center + df_minus) / (h^2);
end

function integrand = computeIntegrand(sw, d2fw_dsw2_func, mob_total_func)
% COMPUTEINTEGRAND Compute integrand: f''(s) / λ_T(s)
    
    % Compute second derivative
    d2fw_dsw2 = d2fw_dsw2_func(sw);
    
    % Compute total mobility
    mob_total = mob_total_func(sw);

    % Compute integrand
    integrand = d2fw_dsw2 ./ mob_total;
    
    % Handle potential NaN or Inf values
    integrand(isnan(integrand) | isinf(integrand)) = 0;
end

function [gamma, details] = computeGammaIntegral(...
    lower_bound, sw_front, integrand_func, n_points, method, verbose)
% COMPUTEGAMMAINTEGRAL Compute gamma integral for given front saturation
    
    % Create saturation grid
    sw_grid = linspace(lower_bound, sw_front, n_points);
    
    % Evaluate functions on grid
    integrand_values = integrand_func(sw_grid);
    
    % Compute integral based on selected method
    switch method
        case 'trapz'
            % Trapezoidal integration
            gamma = trapz(sw_grid, integrand_values);
            
            if verbose
                fprintf('  Trapezoidal integration result: %.6f\n', gamma);
            end
            
        case 'quad'
            % Adaptive quadrature
            gamma = integral(integrand_func, lower_bound, sw_front, 'AbsTol', 1e-10, 'RelTol', 1e-6, 'ArrayValued', true);
            
            if verbose
                fprintf('  Adaptive quadrature result: %.6f\n', gamma);
            end
    end
    
    % Store details
    details = struct();
    details.sw_grid = sw_grid;
    details.integrand_values = integrand_values;
    details.integral_value = gamma;
    details.integration_method = method;
    details.n_points = n_points;
    
    % Compute additional statistics
    details.integrand_max = max(integrand_values);
    details.integrand_min = min(integrand_values);
    details.integrand_mean = mean(integrand_values);
    details.integrand_std = std(integrand_values);  
    
    if verbose
        fprintf('  Integrand statistics:\n');
        fprintf('    Max: %.6e, Min: %.6e, Mean: %.6e\n', ...
                details.integrand_max, details.integrand_min, details.integrand_mean);
    end
end

% =========================================================================
% VISUALIZATION FUNCTIONS (modified for process-specific integration bounds)
% =========================================================================

function plotGammaIntegrandMultiLayer(integration_details, verbose, combine_plots, layer_names)
% PLOTGAMMAINTEGRANDMULTILAYER Create plots for gamma integrand analysis for multiple layers
    
    n_layers = integration_details.n_layers;
    compute_mode = integration_details.Compute;
    
    % Determine figure layout
    if combine_plots && n_layers > 1
        % Combined figure with subplots
        createCombinedGammaFigure(integration_details, verbose, layer_names);
    else
        % Separate figures for each layer
        for layer = 1:n_layers
            createSingleLayerGammaFigure(integration_details, layer, verbose, layer_names);
        end
    end
end

function createCombinedGammaFigure(integration_details, verbose, layer_names)
% CREATECOMBINEDGAMMAFIGURE Create figure with subplots for multiple layers
    
    n_layers = integration_details.n_layers;
    compute_mode = integration_details.Compute;
    process = integration_details.Process;
    lower_bound_label = integration_details.lower_bound_label;
    
    % Determine subplot layout
    n_cols = min(3, n_layers);
    n_rows = ceil(n_layers / n_cols);
    
    % Create figure
    fig = figure('Position', [50, 50, 1600, 900], 'Color', 'white');
    
    % Create title based on compute mode and process
    switch compute_mode
        case 'both'
            title_str = sprintf('Gamma Parameter Analysis (Miscible & Immiscible) - %s Process', process);
        case 'miscible'
            title_str = sprintf('Gamma Parameter Analysis (Miscible Only) - %s Process', process);
        case 'immiscible'
            title_str = sprintf('Gamma Parameter Analysis (Immiscible Only) - %s Process', process);
    end
    sgtitle(title_str, 'FontSize', 16, 'FontWeight', 'bold');
    
    for layer = 1:n_layers
        % Create subplot
        ax = subplot(n_rows, n_cols, layer);
        
        % Plot layer results
        plotLayerGammaResults(ax, integration_details, layer);
        
        % Add layer title
        if ~isempty(layer_names) && length(layer_names) >= layer
            title_str = layer_names{layer};
        else
            title_str = sprintf('Layer %d', layer);
        end
        
        % Add gamma values to title based on compute mode
        switch compute_mode
            case 'both'
                title_str = sprintf('%s\nγ_m = %.4f, γ_i = %.4f', ...
                    title_str, integration_details.gamma_miscible(layer), ...
                    integration_details.gamma_immiscible(layer));
            case 'miscible'
                title_str = sprintf('%s\nγ_m = %.4f', ...
                    title_str, integration_details.gamma_miscible(layer));
            case 'immiscible'
                title_str = sprintf('%s\nγ_i = %.4f', ...
                    title_str, integration_details.gamma_immiscible(layer));
        end
        title(ax, title_str, 'FontSize', 12);
        
        % Add x-axis label with integration bounds
        xlabel(ax, sprintf('Water Saturation (s_w)\nIntegration from %s = %.3f', ...
               lower_bound_label, integration_details.lower_bound(layer)));
        
        % Add grid
        grid(ax, 'on');
        box(ax, 'on');
    end
    
    % Add overall legend to first subplot
    axes(subplot(n_rows, n_cols, 1));
    legend('show', 'Location', 'best');
    
    % Add summary text box
    annotation_str = createGammaSummaryAnnotation(integration_details);
    annotation(fig, 'textbox', [0.02, 0.02, 0.3, 0.1], ...
        'String', annotation_str, 'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 9);
    
    if verbose
        fprintf('  Combined plot generated for %d layers.\n', n_layers);
    end
end

function createSingleLayerGammaFigure(integration_details, layer, verbose, layer_names)
% CREATESINGLELAYERGAMMAFIGURE Create figure for single layer
    
    process = integration_details.Process;
    lower_bound_label = integration_details.lower_bound_label;
    
    % Create figure
    fig = figure('Position', [50, 50, 1400, 800], 'Color', 'white');
    
    % Create axes
    ax = axes('Parent', fig);
    
    % Plot layer results
    plotLayerGammaResults(ax, integration_details, layer);
    
    % Add title
    if ~isempty(layer_names) && length(layer_names) >= layer
        title_str = layer_names{layer};
    else
        title_str = sprintf('Layer %d', layer);
    end
    
    % Add compute mode and process to title
    switch integration_details.Compute
        case 'both'
            title_str = sprintf('%s - Gamma Parameter Analysis (Miscible & Immiscible)', title_str);
        case 'miscible'
            title_str = sprintf('%s - Gamma Parameter Analysis (Miscible)', title_str);
        case 'immiscible'
            title_str = sprintf('%s - Gamma Parameter Analysis (Immiscible)', title_str);
    end
    
    title_str = sprintf('%s - %s Process', title_str, process);
    title(ax, title_str, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add x-axis label with integration bounds
    xlabel(ax, sprintf('Water Saturation (s_w)\nIntegration from %s = %.3f', ...
           lower_bound_label, integration_details.lower_bound(layer)));
    
    % Add results annotation
    annotation_str = createLayerGammaAnnotation(integration_details, layer);
    annotation(fig, 'textbox', [0.02, 0.02, 0.25, 0.15], ...
        'String', annotation_str, 'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    
    if verbose
        fprintf('  Plot generated for layer %d.\n', layer);
    end
end

function plotLayerGammaResults(ax, integration_details, layer)
% PLOTLAYERGAMMARESULTS Plot gamma results for a specific layer
    
    compute_mode = integration_details.Compute;
    process = integration_details.Process;
    lower_bound = integration_details.lower_bound(layer);
    lower_bound_label = integration_details.lower_bound_label;
    
    % Colors
    color_miscible = [0, 0.447, 0.741];  % Blue
    color_immiscible = [0.85, 0.325, 0.098];  % Orange
    
    hold(ax, 'on');
    
    % Plot based on compute mode
    switch compute_mode
        case 'both'
            % Get both sets of data
            details_misc = integration_details.details_miscible{layer};
            details_immisc = integration_details.details_immiscible{layer};
            sw_front_misc = integration_details.sw_front_miscible(layer);
            sw_front_immisc = integration_details.sw_front_immiscible(layer);
            
            % Plot miscible integrand
            plot(ax, details_misc.sw_grid, details_misc.integrand_values, ...
                 '-', 'Color', color_miscible, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Miscible (γ=%.4f)', integration_details.gamma_miscible(layer)));
            
            % Plot immiscible integrand
            plot(ax, details_immisc.sw_grid, details_immisc.integrand_values, ...
                 '-', 'Color', color_immiscible, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Immiscible (γ=%.4f)', integration_details.gamma_immiscible(layer)));
            
            % Mark integration bounds
            xline(ax, lower_bound, 'k--', 'LineWidth', 1, ...
                  'DisplayName', sprintf('%s = %.3f', lower_bound_label, lower_bound));
            xline(ax, sw_front_misc, '--', 'Color', color_miscible, ...
                  'LineWidth', 1, 'DisplayName', 'Miscible Front');
            xline(ax, sw_front_immisc, '--', 'Color', color_immiscible, ...
                  'LineWidth', 1, 'DisplayName', 'Immiscible Front');
            
            % Fill integration areas (from lower_bound to sw_front)
            fill(ax, [lower_bound, details_misc.sw_grid, sw_front_misc], ...
                 [0, details_misc.integrand_values, 0], color_miscible, ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Miscible Area');
            
            fill(ax, [lower_bound, details_immisc.sw_grid, sw_front_immisc], ...
                 [0, details_immisc.integrand_values, 0], color_immiscible, ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Immiscible Area');
            
        case 'miscible'
            % Get miscible data only
            details_misc = integration_details.details_miscible{layer};
            sw_front_misc = integration_details.sw_front_miscible(layer);
            
            % Plot miscible integrand
            plot(ax, details_misc.sw_grid, details_misc.integrand_values, ...
                 '-', 'Color', color_miscible, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Miscible (γ=%.4f)', integration_details.gamma_miscible(layer)));
            
            % Mark integration bounds
            xline(ax, lower_bound, 'k--', 'LineWidth', 1, ...
                  'DisplayName', sprintf('%s = %.3f', lower_bound_label, lower_bound));
            xline(ax, sw_front_misc, '--', 'Color', color_miscible, ...
                  'LineWidth', 1, 'DisplayName', 'Miscible Front');
            
            % Fill integration area (from lower_bound to sw_front)
            fill(ax, [lower_bound, details_misc.sw_grid, sw_front_misc], ...
                 [0, details_misc.integrand_values, 0], color_miscible, ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Miscible Area');
            
        case 'immiscible'
            % Get immiscible data only
            details_immisc = integration_details.details_immiscible{layer};
            sw_front_immisc = integration_details.sw_front_immiscible(layer);
            
            % Plot immiscible integrand
            plot(ax, details_immisc.sw_grid, details_immisc.integrand_values, ...
                 '-', 'Color', color_immiscible, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Immiscible (γ=%.4f)', integration_details.gamma_immiscible(layer)));
            
            % Mark integration bounds
            xline(ax, lower_bound, 'k--', 'LineWidth', 1, ...
                  'DisplayName', sprintf('%s = %.3f', lower_bound_label, lower_bound));
            xline(ax, sw_front_immisc, '--', 'Color', color_immiscible, ...
                  'LineWidth', 1, 'DisplayName', 'Immiscible Front');
            
            % Fill integration area (from lower_bound to sw_front)
            fill(ax, [lower_bound, details_immisc.sw_grid, sw_front_immisc], ...
                 [0, details_immisc.integrand_values, 0], color_immiscible, ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Immiscible Area');
    end
    
    xlabel(ax, 'Water Saturation (s_w)');
    ylabel(ax, 'Integrand: f''''(s) / λ_T(s)');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');
    box(ax, 'on');
    
    % Set appropriate x-axis limits based on process
    if strcmp(process, 'drainage')
        xlim(ax, [0, 1]);
    else
        % For imbibition, show from lower_bound to 1
        xlim(ax, [lower_bound - 0.05, 1]);
    end
    
    hold(ax, 'off');
end

function annotation_str = createGammaSummaryAnnotation(integration_details)
% CREATEGAMMASUMMARYANNOTATION Create summary annotation text for multi-layer
    
    n_layers = integration_details.n_layers;
    compute_mode = integration_details.Compute;
    process = integration_details.Process;
    lower_bound_label = integration_details.lower_bound_label;
    
    switch compute_mode
        case 'both'
            annotation_str = sprintf('Multi-Layer Gamma Summary (Both):\n');
            annotation_str = sprintf('%sProcess: %s\n', annotation_str, process);
            annotation_str = sprintf('%sIntegration from: %s\n\n', annotation_str, lower_bound_label);
            
            for layer = 1:min(5, n_layers)
                annotation_str = sprintf('%sL%d: γ_m=%.4f, γ_i=%.4f\n', ...
                    annotation_str, layer, integration_details.gamma_miscible(layer), ...
                    integration_details.gamma_immiscible(layer));
            end
            
            if n_layers > 5
                annotation_str = sprintf('%s... and %d more layers', annotation_str, n_layers - 5);
            end
            
            % Add statistics
            avg_gamma_m = mean(integration_details.gamma_miscible);
            avg_gamma_i = mean(integration_details.gamma_immiscible);
            avg_ratio = mean(integration_details.gamma_miscible ./ integration_details.gamma_immiscible);
            
            annotation_str = sprintf('%s\nAverages:\nγ_m: %.4f\nγ_i: %.4f\nRatio: %.4f', ...
                annotation_str, avg_gamma_m, avg_gamma_i, avg_ratio);
            
        case 'miscible'
            annotation_str = sprintf('Multi-Layer Gamma Summary (Miscible):\n');
            annotation_str = sprintf('%sProcess: %s\n', annotation_str, process);
            annotation_str = sprintf('%sIntegration from: %s\n\n', annotation_str, lower_bound_label);
            
            for layer = 1:min(5, n_layers)
                annotation_str = sprintf('%sL%d: γ=%.4f\n', ...
                    annotation_str, layer, integration_details.gamma_miscible(layer));
            end
            
            if n_layers > 5
                annotation_str = sprintf('%s... and %d more layers', annotation_str, n_layers - 5);
            end
            
            % Add statistics
            avg_gamma = mean(integration_details.gamma_miscible);
            min_gamma = min(integration_details.gamma_miscible);
            max_gamma = max(integration_details.gamma_miscible);
            
            annotation_str = sprintf('%s\nStatistics:\nAverage: %.4f\nMin: %.4f\nMax: %.4f', ...
                annotation_str, avg_gamma, min_gamma, max_gamma);
            
        case 'immiscible'
            annotation_str = sprintf('Multi-Layer Gamma Summary (Immiscible):\n');
            annotation_str = sprintf('%sProcess: %s\n', annotation_str, process);
            annotation_str = sprintf('%sIntegration from: %s\n\n', annotation_str, lower_bound_label);
            
            for layer = 1:min(5, n_layers)
                annotation_str = sprintf('%sL%d: γ=%.4f\n', ...
                    annotation_str, layer, integration_details.gamma_immiscible(layer));
            end
            
            if n_layers > 5
                annotation_str = sprintf('%s... and %d more layers', annotation_str, n_layers - 5);
            end
            
            % Add statistics
            avg_gamma = mean(integration_details.gamma_immiscible);
            min_gamma = min(integration_details.gamma_immiscible);
            max_gamma = max(integration_details.gamma_immiscible);
            
            annotation_str = sprintf('%s\nStatistics:\nAverage: %.4f\nMin: %.4f\nMax: %.4f', ...
                annotation_str, avg_gamma, min_gamma, max_gamma);
    end
end

function annotation_str = createLayerGammaAnnotation(integration_details, layer)
% CREATELAYERGAMMAANNOTATION Create annotation for single layer
    
    compute_mode = integration_details.Compute;
    process = integration_details.Process;
    lower_bound = integration_details.lower_bound(layer);
    lower_bound_label = integration_details.lower_bound_label;
    
    switch compute_mode
        case 'both'
            annotation_str = sprintf(['Layer %d Results (Both):\n\n' ...
                                     'Process: %s\n' ...
                                     'Integration from: %s = %.4f\n' ...
                                     'Miscible Front: %.4f\n' ...
                                     'Immiscible Front: %.4f\n\n' ...
                                     'Gamma Parameters:\n' ...
                                     '  γ_miscible: %.6f\n' ...
                                     '  γ_immiscible: %.6f\n' ...
                                     '  Ratio (M/I): %.4f\n\n' ...
                                     'Integrand Statistics:\n' ...
                                     '  Max (M): %.2e\n' ...
                                     '  Mean (M): %.2e'], ...
                                    layer, ...
                                    process, ...
                                    lower_bound_label, lower_bound, ...
                                    integration_details.sw_front_miscible(layer), ...
                                    integration_details.sw_front_immiscible(layer), ...
                                    integration_details.gamma_miscible(layer), ...
                                    integration_details.gamma_immiscible(layer), ...
                                    integration_details.gamma_miscible(layer) / integration_details.gamma_immiscible(layer), ...
                                    integration_details.details_miscible{layer}.integrand_max, ...
                                    integration_details.details_miscible{layer}.integrand_mean);
                                
        case 'miscible'
            annotation_str = sprintf(['Layer %d Results (Miscible):\n\n' ...
                                     'Process: %s\n' ...
                                     'Integration from: %s = %.4f\n' ...
                                     'Miscible Front: %.4f\n\n' ...
                                     'Gamma Parameter:\n' ...
                                     '  γ_miscible: %.6f\n\n' ...
                                     'Integrand Statistics:\n' ...
                                     '  Max: %.2e\n' ...
                                     '  Mean: %.2e'], ...
                                    layer, ...
                                    process, ...
                                    lower_bound_label, lower_bound, ...
                                    integration_details.sw_front_miscible(layer), ...
                                    integration_details.gamma_miscible(layer), ...
                                    integration_details.details_miscible{layer}.integrand_max, ...
                                    integration_details.details_miscible{layer}.integrand_mean);
                                
        case 'immiscible'
            annotation_str = sprintf(['Layer %d Results (Immiscible):\n\n' ...
                                     'Process: %s\n' ...
                                     'Integration from: %s = %.4f\n' ...
                                     'Immiscible Front: %.4f\n\n' ...
                                     'Gamma Parameter:\n' ...
                                     '  γ_immiscible: %.6f\n\n' ...
                                     'Integrand Statistics:\n' ...
                                     '  Max: %.2e\n' ...
                                     '  Mean: %.2e'], ...
                                    layer, ...
                                    process, ...
                                    lower_bound_label, lower_bound, ...
                                    integration_details.sw_front_immiscible(layer), ...
                                    integration_details.gamma_immiscible(layer), ...
                                    integration_details.details_immiscible{layer}.integrand_max, ...
                                    integration_details.details_immiscible{layer}.integrand_mean);
    end
end