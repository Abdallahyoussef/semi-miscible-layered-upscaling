function [double_int_values, integration_details] = ...
    computeDoubleIntegralTerm(Swr, Snwr, sw_front, sw_lower_array, sw_upper_array, ...
                              fw_func, dfw_dsw_func, mob_t_func, varargin)
% COMPUTEDOUBLEINTEGRALTERM2 Compute double integral term for second stage
%
% DESCRIPTION:
%   Vectorized version that computes the double integral term for multiple layers:
%   ∫_{sw_front}^{sw} [ f''(s) / f'(s)^3 * ∫_{lower_bound}^{s} [ f''(ξ) / λ_T(ξ) ] dξ ] ds
%
%   For imbibition: lower_bound = 1-Snwr
%   For drainage: lower_bound = Swr
%
% SYNTAX:
%   For single layer:
%   double_int_value = computeDoubleIntegralTerm2(Swr, Snwr, sw_lower, sw_upper, ...
%                                                 fw_func, dfw_dsw_func, mob_t_func)
%
%   For multiple layers:
%   double_int_values = computeDoubleIntegralTerm2(Swr, Snwr, sw_lower_array, sw_upper_array, ...
%                                                  fw_func, dfw_dsw_func, mob_t_func)
%
%   [double_int_values, integration_details] = computeDoubleIntegralTerm2(...)
%
% INPUTS:
%   Swr             - Residual water saturation (scalar or vector)
%   Snwr            - Residual non-wetting saturation (scalar or vector)
%   sw_front   - Flood front saturation (scalar or vector)
%   sw_lower_array  - Lower bounds for outer integral (scalar or vector)
%   sw_upper_array  - Upper bounds for outer integral (scalar or vector)
%   fw_func         - Function handle or cell array: fw(sw)
%   dfw_dsw_func    - Function handle or cell array: dfw_dsw(sw)
%   mob_t_func      - Function handle or cell array: λ_T(sw)
%
% OPTIONAL INPUTS:
%   layer_idx    - (Optional) Layer index for multi-layer systems (default: 1)
%
% OPTIONAL PARAMETERS:
%   'Process'           - 'drainage' or 'imbibition' (default: 'imbibition')
%   'InnerPoints'       - Number of points for inner integral (default: 200)
%   'OuterPoints'       - Number of points for outer integral (default: 200)
%   'IntegrationMethod' - 'trapz' or 'integral' (default: 'integral')
%   'SecondDerivMethod' - 'complex' or 'finite' (default: 'finite')
%   'ComplexStepSize'   - Step size for complex step derivative (default: 1e-20)
%   'FiniteDiffStep'    - Step size for finite difference (default: 1e-6)
%   'Verbose'           - Display progress (default: true)
%   'Precompute'        - Precompute inner integrals (default: true)
%
% OUTPUTS:
%   double_int_values   - Double integral values (scalar or vector)
%   integration_details - Structure with detailed results
%

% =========================================================================
% PARSE INPUTS
% =========================================================================

% Default parameters
defaults.Process = 'imbibition';
defaults.InnerPoints = 200;
defaults.OuterPoints = 200;
defaults.IntegrationMethod = 'integral';
defaults.SecondDerivMethod = 'finite';
defaults.ComplexStepSize = 1e-20;
defaults.FiniteDiffStep = 1e-6;
defaults.Verbose = true;
defaults.Precompute = true;

% Check if layer_idx is provided as numeric argument
layer_idx = [];
if ~isempty(varargin) && isscalar(varargin{1}) && isPositiveIntegerArray(varargin{1})
    layer_idx = varargin{1};
    varargin = varargin(2:end);
end

% Parse name-value pairs
p = inputParser;
addRequired(p, 'Swr', @(x) isnumeric(x));
addRequired(p, 'Snwr', @(x) isnumeric(x));
addRequired(p, 'sw_lower_array', @(x) isnumeric(x));
addRequired(p, 'sw_upper_array', @(x) isnumeric(x));
addRequired(p, 'fw_func');
addRequired(p, 'dfw_dsw_func');
addRequired(p, 'mob_t_func');

fields = fieldnames(defaults);
for i = 1:length(fields)
    addParameter(p, fields{i}, defaults.(fields{i}));
end

parse(p, Swr, Snwr, sw_lower_array, sw_upper_array, ...
      fw_func, dfw_dsw_func, mob_t_func, varargin{:});
params = p.Results;

% Validate Process parameter
valid_process = {'drainage', 'imbibition', 'Drainage', 'Imbibition'};
params.Process = validatestring(params.Process, valid_process);

if strcmp(params.Process, 'Drainage') || strcmp(params.Process, 'drainage')
    params.Process = 'drainage';
else
    params.Process = 'imbibition';
end

% Validate integration method
valid_methods = {'trapz', 'integral'};
params.IntegrationMethod = validatestring(params.IntegrationMethod, valid_methods);

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
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of fractional flow functions must match derivative functions');

        % Validate that mobility functions match
        elseif ~iscell(mob_t_func) || length(mob_t_func) ~= n_layers
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of mobility functions must match fractional flow functions');

        % Validate that residual water saturations match
        elseif ~isvector(Swr) || length(Swr) ~= n_layers
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of residual water saturations must match number of layers');

        % Validate that residual non-wetting saturations match
        elseif ~isvector(Snwr) || length(Snwr) ~= n_layers
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of residual non-wetting saturations must match number of layers');

        end
    
    % Check if calculations are required for specific layers
    elseif ~isPositiveIntegerArray(layer_idx) || max(layer_idx) > n_layers
        error('computeDoubleIntegralTerm:InvalidInput', ...
            'layer_idx must be positive integer and less than number of layers')

    else   % only specific layers are required
        n_layers = length(layer_idx);

        % Validate that residual saturations match
        if ~isvector(Swr) || length(Swr) ~= n_layers
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of residual water saturations must match layer_idx');
        end
        
        if ~isvector(Snwr) || length(Snwr) ~= n_layers
            error('computeDoubleIntegralTerm:InvalidInput', ...
                  'Number of residual non-wetting saturations must match layer_idx');
        end

        % Extract specific layer functions
        fw_func = fw_func(layer_idx);
        dfw_dsw_func = dfw_dsw_func(layer_idx);
        mob_t_func = mob_t_func(layer_idx);
    end

    % Validate saturation arrays
    if ~isvector(sw_lower_array) || length(sw_lower_array) ~= n_layers
        error('computeDoubleIntegralTerm:InvalidInput', ...
              'Number of lower bounds must match number of layers');
    end
    
    if ~isvector(sw_upper_array) || length(sw_upper_array) ~= n_layers
        error('computeDoubleIntegralTerm:InvalidInput', ...
              'Number of upper bounds must match number of layers');
    end

else
    n_layers = 1;
    
    % Convert to cell arrays for uniform handling
    fw_func = {fw_func};
    dfw_dsw_func = {dfw_dsw_func};
    mob_t_func = {mob_t_func};
    
end

% Ensure saturations are vectors
if isscalar(Swr)
    Swr = Swr(:);
end
if isscalar(Snwr)
    Snwr = Snwr(:);
end
if isscalar(sw_lower_array)
    sw_lower_array = sw_lower_array(:);
end
if isscalar(sw_upper_array)
    sw_upper_array = sw_upper_array(:);
end


% Check if input endpoint saturations are valid
if min(Swr) < 0 || max(Swr) >= 1 || min(Snwr) < 0 || max(Snwr) >= 1 || max(Swr + Snwr) >= 1
        error('computeDoubleIntegralTerm:InvalidResiduals', ...
              ['Invalid residual saturations. Input residual saturations must be between zero and 1, ' ...
              '& their sum less than or equal 1']);
end

% Check if front saturations are valid
if any((sw_front - Swr) < 0) || any((sw_front - (1 - Snwr)) > 0)
        error('computeDoubleIntegralTerm:InvalidFrontSaturation', ...
              ['Invalid front saturations. Input front saturations must be between Swr and 1-Snwr, ' ...
              '& their sum less than or equal 1']);
end

% Calculate lower bound for inner integral based on process
if strcmp(params.Process, 'drainage')
    inner_lower = Swr;  % For drainage: integrate from Swr
    inner_lower_label = 'Swr';
else
    inner_lower = 1 - Snwr;  % For imbibition: integrate from 1-Snwr
    inner_lower_label = '1-Snwr';
end

% Validate saturation bounds for all layers
parfor layer = 1:n_layers
    
    % For imbibition, check that saturations are less than inner_lower
    if strcmp(params.Process, 'imbibition')
        if sw_lower_array(layer) > inner_lower(layer) || ...
           sw_upper_array(layer) > inner_lower(layer)  || ...
           sw_lower_array(layer) < sw_front(layer) || ...
           sw_upper_array(layer) < sw_front(layer)
            error('computeDoubleIntegralTerm:InvalidBoundsImbibition', ...
                  'Layer %d: For imbibition, saturations must be less than or equal %s = %.4f & more than or equal sw_front = %.4f', ...
                  layer, inner_lower_label, inner_lower(layer), sw_front(layer));
        end
    % For drainage, check that saturations are more than inner_lower
    else
        if sw_lower_array(layer) < inner_lower(layer) || ...
           sw_upper_array(layer) < inner_lower(layer) || ...
           sw_lower_array(layer) > sw_front(layer) || ...
           sw_upper_array(layer) > sw_front(layer)
            error('computeDoubleIntegralTerm:InvalidBoundsImbibition', ...
                  'Layer %d: For drainage, saturations must be more than or equal %s = %.4f & less than or equal sw_front = %.4f', ...
                  layer, inner_lower_label, inner_lower(layer), sw_front(layer));
        end
    end
end

if params.Verbose
    fprintf('\n=== Computing Double Integral Term ===\n');
    fprintf('Number of layers: %d\n', n_layers);
    fprintf('Process: %s\n', params.Process);
    fprintf('Integration from inner bound: %s\n', inner_lower_label);
    fprintf('Integration method: %s\n', params.IntegrationMethod);
    if n_layers > 1
        fprintf('Inner bounds (%s): %s\n', inner_lower_label, mat2str(inner_lower', 3));
    end
end

% =========================================================================
% PROCESS EACH LAYER
% =========================================================================

% Initialize outputs
double_int_values = zeros(n_layers, 1);

% Prepare results structure
integration_details = struct();
integration_details.n_layers = n_layers;
integration_details.Process = params.Process;
integration_details.Swr = Swr;
integration_details.Snwr = Snwr;
integration_details.inner_lower = inner_lower;
integration_details.inner_lower_label = inner_lower_label;
integration_details.sw_lower_array = sw_lower_array;
integration_details.sw_upper_array = sw_upper_array;
integration_details.double_int_values = zeros(n_layers, 1);
integration_details.details = cell(n_layers, 1);
integration_details.methods = struct(...
    'integration_method', params.IntegrationMethod, ...
    'second_derivative_method', params.SecondDerivMethod);

% Looping on layers to compute double integrals
for layer = 1:n_layers
    if params.Verbose
        fprintf('\n--- Processing Layer %d ---\n', layer);
        fprintf('Inner bound (%s): %.4f\n', inner_lower_label, inner_lower(layer));
        fprintf('Outer integral: [%.4f, %.4f]\n', ...
                sw_lower_array(layer), sw_upper_array(layer));
    end
    
    if sw_lower_array(layer) ~= sw_upper_array(layer)
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
        
        % Create function for second derivative
        d2fw_dsw2_func = createSecondDerivativeFunction(...
            current_fw_func, params.SecondDerivMethod, ...
            params.ComplexStepSize, params.FiniteDiffStep, params.Verbose);
        
        % =====================================================================
        % CREATE INNER INTEGRAND FUNCTION
        % =====================================================================
        
        if params.Verbose
            fprintf('Creating inner integrand: f''''(ξ) / λ_T(ξ)\n');
        end
        
        inner_integrand_func = @(xi) safeDivision(d2fw_dsw2_func(xi), current_mob_t_func(xi));
        
        % =====================================================================
        % CREATE INNER INTEGRAL FUNCTION G(s)
        % =====================================================================
    
           if params.Precompute
                % Precompute G(s) on a grid for efficiency
                s_grid = linspace(sw_lower_array(layer), sw_upper_array(layer), params.OuterPoints);
                G_values = zeros(size(s_grid));
                
                parfor i = 1:length(s_grid)
                    s = s_grid(i);
                    if s ~= inner_lower(layer)
                        G_values(i) = computeIntegral(inner_integrand_func, ...
                                                     inner_lower(layer), s, ...
                                                     params.InnerPoints, params.IntegrationMethod);
                    else
                        G_values(i) = 0;
                    end
                end
                
                % Create interpolation function for G(s)
                G_func = @(s) interp1(s_grid, G_values, s, 'pchip', 'extrap');
                
                if params.Verbose
                    fprintf('  G(s) precomputed at %d points\n', length(s_grid));
                end
            else
                % Create G function that computes on the fly
                G_func = @(s) computeInnerIntegral(s, inner_lower(layer), ...
                                                  inner_integrand_func, ...
                                                  params.InnerPoints, params.IntegrationMethod);
            end
            
            % =====================================================================
            % CREATE OUTER INTEGRAND FUNCTION
            % =====================================================================
            
            if params.Verbose
                fprintf('Creating outer integrand: f''''(s) / f''(s)^3 * G(s)\n');
            end
            
            outer_integrand_func = @(s) safeDivision(d2fw_dsw2_func(s), ...
                                                   current_dfw_func(s).^3) .* G_func(s);
            
            % =====================================================================
            % COMPUTE DOUBLE INTEGRAL
            % =====================================================================
            
            if params.Verbose
                fprintf('Computing double integral...\n');
            end
            
            double_int_layer = computeIntegral(outer_integrand_func, ...
                                              sw_lower_array(layer), sw_upper_array(layer), ...
                                              params.OuterPoints, params.IntegrationMethod);
    
        % Store results
        double_int_values(layer) = double_int_layer;
        integration_details.double_int_values(layer) = double_int_layer;
        
        % Store layer details
        layer_details = struct();
        layer_details.double_int_value = double_int_layer;
        if params.Precompute
            layer_details.s_grid = s_grid;
            layer_details.G_values = G_values;
        end
        integration_details.details{layer} = layer_details;
    else
        double_int_layer = 0;
        double_int_values(layer) = double_int_layer;
    end
    
    if params.Verbose
        fprintf('  Double integral value: %.6e\n', double_int_layer);
    end
end

% =========================================================================
% DISPLAY FINAL RESULTS
% =========================================================================

if params.Verbose
    fprintf('\n=== Double Integral Results ===\n');
    fprintf('Process: %s\n', params.Process);
    
    if n_layers == 1
        fprintf('Double integral value: %.6e\n', double_int_values(1));
    else
        fprintf('Layer\tDouble Integral Value\n');
        fprintf('-----\t---------------------\n');
        for layer = 1:n_layers
            fprintf('%d\t%.6e\n', layer, double_int_values(layer));
        end
        
        fprintf('\nSummary Statistics:\n');
        fprintf('  Average: %.6e\n', mean(double_int_values));
        fprintf('  Min: %.6e\n', min(double_int_values));
        fprintf('  Max: %.6e\n', max(double_int_values));
        fprintf('  Std: %.6e\n', std(double_int_values));
    end
end

if params.Verbose
    fprintf('\n=== Double Integral Computation Complete ===\n\n');
end

end % Main function

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function d2fw_dsw2_func = createSecondDerivativeFunction(...
    fw_func, method, complex_step, finite_step, verbose)
% CREATESECONDDERIVATIVEFUNCTION Create function for second derivative
    
    switch method
        case 'complex'
            % Complex step method
            if verbose
                fprintf('  Using complex step method (step = %.1e)\n', complex_step);
            end
            
            d2fw_dsw2_func = @(sw) imag(fw_func(sw + complex_step*1i)) / complex_step^2;
            
        case 'finite'
            % Finite difference method
            if verbose
                fprintf('  Using finite difference method (step = %.1e)\n', finite_step);
            end
            
            d2fw_dsw2_func = @(sw) (fw_func(sw + finite_step) - 2*fw_func(sw) + ...
                                   fw_func(sw - finite_step)) / finite_step^2;
            
        otherwise
            error('computeDoubleIntegralTerm2:InvalidMethod', ...
                  'Unknown second derivative method: %s', method);
    end
end

function integral_value = computeIntegral(integrand_func, lower, upper, ...
                                         n_points, method)
% COMPUTEINTEGRAL Compute integral using specified method
    
    switch method
        case 'trapz'
            % Trapezoidal integration
            s_grid = linspace(lower, upper, n_points);
            integrand_vals = integrand_func(s_grid);
            integral_value = trapz(s_grid, integrand_vals);
            
        case 'integral'
            % Adaptive quadrature
            integral_value = integral(integrand_func, lower, upper, ...
                                     'AbsTol', 1e-10, 'RelTol', 1e-6);
    end
end

function G_value = computeInnerIntegral(s, inner_lower, inner_integrand_func, ...
                                       n_points, method)
% COMPUTEINNERINTEGRAL Compute inner integral G(s) on the fly

    % Initialize output
    G_value = zeros(size(s));

    % Process each s value
    parfor i = 1:numel(s)
        s_input = s(i);
    
        if s_input <= inner_lower
            G_value(i) = 0;
        else
            G_value(i) = computeIntegral(inner_integrand_func, inner_lower, s_input, ...
                         n_points, method);
        end
    end
    
    % Reshape to match input
    G_value = reshape(G_value, size(s));
    

end

function result = safeDivision(numerator, denominator)
% SAFEDIVISION Safe division handling zeros and singularities
    
    result = numerator ./ denominator;
    
    % Handle division by zero or near-zero
    mask = abs(denominator) < eps;
    if any(mask)
        % Use regularization for near-zero denominators
        denominator_reg = denominator;
        denominator_reg(mask) = sign(denominator(mask)) * eps;
        result = numerator ./ denominator_reg;
    end
    
    % % Handle NaN and Inf
    % result(isnan(result) | isinf(result)) = 0;
end