function varargout = relativePermeability(model_type, varargin)
% RELATIVEPERMEABILITY Create relative permeability and related functions
%
% DESCRIPTION:
%   This function creates function handles for relative permeability and
%   related properties for two-phase flow in porous media. Supports multiple
%   layers with different relative permeability models.
%
% SYNTAX:
%   Basic usage (single layer):
%   [krw_func, krnw_func] = relativePermeability(model_type, parameters)
%
%   With viscosity inputs (extended outputs):
%   [krw_func, krnw_func, mob_w_func, mob_nw_func, mob_t_func, ...
%    fw_func, dfw_dsw_func] = relativePermeability(model_type, parameters, mu_w, mu_nw)
%
%   Multiple layers:
%   func_handles = relativePermeability(model_type, parameters, ...)
%   where model_type and parameters can be cell arrays for multiple layers
%
% INPUTS:
%   model_type - String or cell array of strings specifying correlation models
%       'Brooks-Corey' or 'BC' - Brooks-Corey model
%       'Van-Genuchten' or 'VG' - Van Genuchten model
%       'LET' - LET correlation model
%       'Chierici'- Chierici model
%
%   parameters - Model parameters (structure, cell array of structures, or individual params)
%
%   mu_w - (Optional) Viscosity of water [Pa·s]
%   mu_nw - (Optional) Viscosity of CO2 [Pa·s]
%
% OUTPUTS:
%   krw_func   - Function handle: krw(sw) or krw(sw, layer_idx)
%   krnw_func  - Function handle: krnw(sw) or krnw(sw, layer_idx)
%   mob_w_func - (Optional) Water mobility function
%   mob_nw_func - (Optional) CO2 mobility function
%   mob_t_func  - (Optional) Total mobility function
%   fw_func     - (Optional) Fractional flow function
%   dfw_dsw_func - (Optional) Fractional flow derivative function
%
%   For multiple layers: function handles accept layer index as second argument
%
% EXAMPLES:
%   % Example 1: Single layer with extended outputs
%   mu_w = 1e-3; mu_nw = 6e-5;
%   params.lambda = 2.0; params.Swr = 0.2; params.Snwr = 0.15;
%   params.Krw_end = 0.3;  params.Krnw_end = 0.9;
%   [krw_func, krnw_func, mob_w_func, mob_nw_func, mob_t_func, fw_func, dfw_dsw_func] = ...
%       relativePermeability('BC', params, mu_w, mu_nw);
%
%   % Example 2: Multiple layers
%   models = {'BC', 'VG', 'LET', 'Chierici'};
%   params{1} = struct('lambda', 2.0, 'Swr', 0.2, 'Snwr', 0.15, 'Krw_end', 0.3, 'Krnw_end', 0.9);
%   params{2} = struct('m', 0.5, 'n', 2.0, 'Swr', 0.25, 'Snwr', 0.10, 'Krw_end', 0.3, 'Krnw_end', 0.9);
%   params{3} = struct('Lw',1.0,'Ew',1.5,'Tw',2.0,'Lnw',1.2,'Enw',1.0,'Tnw',1.8,'Swr',0.15,'Snwr',0.20, 'Krw_end', 0.3, 'Krnw_end', 0.9);
%   params{4} = struct('A',1.0,'L',1.5,'B',2.0,'M',1.2,'Swr',0.15,'Snwr',0.20, 'Krw_end', 0.27, 'Krnw_end', 0.95);
%   funcs = relativePermeability(models, params, mu_w, mu_nw);
%
% AUTHOR: AbdAllah AbdRabElNaby Youssef ElSaieed
% DATE: 28 December 2025

% =========================================================================
% INPUT VALIDATION AND PARSING
% =========================================================================

% Check minimum number of arguments
narginchk(2, 4);

% Initialize flags
multi_layer = false;
viscosities_provided = false;
mu_w = [];
mu_nw = [];

% Check if viscosities are provided
if nargin >= 4
    viscosities_provided = true;
    mu_w = varargin{end-1};
    mu_nw = varargin{end};
    
    % Validate viscosities
    if mu_w <= 0 || mu_nw <= 0
        error('relativePermeability:InvalidViscosity', ...
              'Viscosities must be positive');
    end
    
    % Adjust varargin to remove viscosities from parameters
    varargin = varargin(1:end-2);
end

% Check for multiple layers
if iscell(model_type) || (isstruct(varargin{1}) && isfield(varargin{1}, 'models'))
    multi_layer = true;
    n_layers = parseMultiLayerInput(model_type, varargin{1});
else
    n_layers = 1;
end

% Parse parameters based on layer configuration
if multi_layer
    [model_types, param_structs] = parseMultiLayerParameters(model_type, varargin{1}, n_layers);
else
    model_types = {model_type};
    param_structs = {parseSingleLayerParameters(model_type, varargin{:})};
end

% =========================================================================
% CREATE FUNCTION HANDLES
% =========================================================================

if multi_layer
    % Create function handles for multiple layers
    [krw_func, krnw_func] = createMultiLayerFunctions(model_types, param_structs);
    
    if viscosities_provided
        [mob_w_func, mob_nw_func, mob_t_func, fw_func, dfw_dsw_func] = ...
            createMultiLayerExtendedFunctions(model_types, param_structs, mu_w, mu_nw);
    end
else
    % Single layer functions
    [krw_func, krnw_func] = createSingleLayerFunctions(model_types{1}, param_structs{1});
    
    if viscosities_provided
        [mob_w_func, mob_nw_func, mob_t_func, fw_func, dfw_dsw_func] = ...
            createSingleLayerExtendedFunctions(model_types{1}, param_structs{1}, mu_w, mu_nw);
    end
end

% =========================================================================
% SET OUTPUTS
% =========================================================================

if nargout == 0
    % Display information if no outputs requested
    displayModelInfo(model_types, param_structs, viscosities_provided);
    
elseif nargout == 1 && multi_layer
    % Package all functions in a structure for multiple layers
    varargout{1} = struct(...
        'krw_func', krw_func, ...
        'krnw_func', krnw_func);
    
    if viscosities_provided
        varargout{1}.mob_w_func = mob_w_func;
        varargout{1}.mob_nw_func = mob_nw_func;
        varargout{1}.mob_t_func = mob_t_func;
        varargout{1}.fw_func = fw_func;
        varargout{1}.dfw_dsw_func = dfw_dsw_func;
    end
    
else
    % Standard output assignment
    varargout{1} = krw_func;
    varargout{2} = krnw_func;
    
    if nargout >= 3 && viscosities_provided
        varargout{3} = mob_w_func;
    end
    if nargout >= 4 && viscosities_provided
        varargout{4} = mob_nw_func;
    end
    if nargout >= 5 && viscosities_provided
        varargout{5} = mob_t_func;
    end
    if nargout >= 6 && viscosities_provided
        varargout{6} = fw_func;
    end
    if nargout >= 7 && viscosities_provided
        varargout{7} = dfw_dsw_func;
    end
end

end % Main function

% =========================================================================
% HELPER FUNCTIONS - INPUT PARSING
% =========================================================================

function n_layers = parseMultiLayerInput(model_type, parameters)
% PARSEMULTILAYERINPUT Determine number of layers from input
    
    if iscell(model_type)
        n_layers = length(model_type);
    elseif isstruct(parameters) && isfield(parameters, 'models')
        n_layers = length(parameters.models);
    else
        error('relativePermeability:InvalidMultiLayerInput', ...
              'For multiple layers, provide cell arrays or structure with ''models'' field');
    end
end

function [model_types, param_structs] = parseMultiLayerParameters(model_type, parameters, n_layers)
% PARSEMULTILAYERPARAMETERS Parse parameters for multiple layers
    
    model_types = cell(1, n_layers);
    param_structs = cell(1, n_layers);
    
    if iscell(model_type)
        % Cell array input format
        model_types = model_type;
        
        if iscell(parameters)
            % Parameters as cell array of structures
            for i = 1:n_layers
                if isstruct(parameters{i})
                    param_structs{i} = parameters{i};
                else
                    param_structs{i} = parseSingleLayerParameters(...
                        model_types{i}, parameters{i});
                end
            end
        elseif isstruct(parameters) && isfield(parameters, 'layer_params')
            % Parameters in structured format
            for i = 1:n_layers
                param_structs{i} = parameters.layer_params{i};
            end
        else
            error('relativePermeability:InvalidMultiLayerParams', ...
                  'Invalid parameter format for multiple layers');
        end
        
    elseif isstruct(parameters) && isfield(parameters, 'models')
        % Structured input format
        for i = 1:n_layers
            model_types{i} = parameters.models{i};
            param_structs{i} = parameters.params{i};
        end
    end
end

function params = parseSingleLayerParameters(model_type, varargin)
% PARSESINGLELAYERPARAMETERS Parse parameters for single layer
    
    if isstruct(varargin{1})
        params = varargin{1};
    else
        % Parse individual parameters based on model type
        model_type = validatestring(upper(model_type), ...
            {'BROOKS-COREY', 'BC', 'VAN-GENUCHTEN', 'VG', 'LET', 'Chierici'});
        
        params = struct();
        switch model_type
            case {'BC', 'BROOKS-COREY'}
                if length(varargin) >= 5
                    params.lambda = varargin{1};
                    params.Swr = varargin{2};
                    params.Snwr = varargin{3};
                    params.Krw_end = varargin{4};
                    params.Krnw_end = varargin{5};
                end
            case {'VG', 'VAN-GENUCHTEN'}
                if length(varargin) >= 6
                    params.m = varargin{1};
                    params.n = varargin{2};
                    params.Swr = varargin{3};
                    params.Snwr = varargin{4};
                    params.Krw_end = varargin{5};
                    params.Krnw_end = varargin{6};
                end
            case 'LET'
                if length(varargin) >= 10
                    params.Lw = varargin{1};
                    params.Ew = varargin{2};
                    params.Tw = varargin{3};
                    params.Lnw = varargin{4};
                    params.Enw = varargin{5};
                    params.Tnw = varargin{6};
                    params.Swr = varargin{7};
                    params.Snwr = varargin{8};
                    params.Krw_end = varargin{9};
                    params.Krnw_end = varargin{10};
                end
            case {'Chierici'}
                if length(varargin) >= 8
                    params.A = varargin{1};
                    params.L = varargin{2};
                    params.B = varargin{3};
                    params.M = varargin{4};
                    params.Swr = varargin{5};
                    params.Snwr = varargin{6};
                    params.Krw_end = varargin{7};
                    params.Krnw_end = varargin{8};
                end
        end
    end
    
    % Validate required parameters
    validateParameters(params, model_type);
end

% =========================================================================
% HELPER FUNCTIONS - FUNCTION CREATION
% =========================================================================

function [krw_func, krnw_func] = createSingleLayerFunctions(model_type, params)
% CREATESINGLELAYERFUNCTIONS Create function handles for single layer
    
    % Standardize model type
    model_type = standardizeModelType(model_type);
    
    % Create function handles
    switch model_type
        case 'BROOKS-COREY'
            krw_func = @(sw) brooksCoreyKrw(sw, params.lambda, params.Swr, params.Snwr, params.Krw_end);
            krnw_func = @(sw) brooksCoreyKrnw(sw, params.lambda, params.Swr, params.Snwr, params.Krnw_end);
            
        case 'VAN-GENUCHTEN'
            krw_func = @(sw) vanGenuchtenKrw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krw_end);
            krnw_func = @(sw) vanGenuchtenKrnw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krnw_end);
            
        case 'LET'
            krw_func = @(sw) letKrw(sw, params.Lw, params.Ew, params.Tw, params.Swr, params.Snwr, params.Krw_end);
            krnw_func = @(sw) letKrnw(sw, params.Lnw, params.Enw, params.Tnw, params.Swr, params.Snwr, params.Krnw_end);

        case 'Chierici'
            krw_func = @(sw) chiericiKrw(sw, params.A, params.L, params.Swr, params.Snwr, params.Krw_end);
            krnw_func = @(sw) chiericiKrnw(sw, params.B, params.M, params.Swr, params.Snwr, params.Krnw_end);
    end
end

function [krw_func, krnw_func] = createMultiLayerFunctions(model_types, param_structs)
% CREATEMULTILAYERFUNCTIONS Create function handles for multiple layers
    
    n_layers = length(model_types);
    
    krw_func = @(sw, layer_idx) evaluateLayerKrw(sw, layer_idx, model_types, param_structs, n_layers);
    krnw_func = @(sw, layer_idx) evaluateLayerKrnw(sw, layer_idx, model_types, param_structs, n_layers);
end

function [mob_w_func, mob_nw_func, mob_t_func, fw_func, dfw_dsw_func] = ...
    createSingleLayerExtendedFunctions(model_type, params, mu_w, mu_nw)
% CREATESINGLELAYEREXTENDEDFUNCTIONS Create extended function handles for single layer
    
    % Get basic relative permeability functions
    [krw_base_func, krnw_base_func] = createSingleLayerFunctions(model_type, params);
    
    % Create mobility functions
    mob_w_func = @(sw) krw_base_func(sw) / mu_w;
    mob_nw_func = @(sw) krnw_base_func(sw) / mu_nw;
    mob_t_func = @(sw) mob_w_func(sw) + mob_nw_func(sw);
    
    % Create fractional flow function
    fw_func = @(sw) mob_w_func(sw) ./ (mob_w_func(sw) + mob_nw_func(sw) + eps);
    
    % Create fractional flow derivative function using complex step differentiation
    dfw_dsw_func = @(sw) complexStepDerivative(fw_func, sw);
end

function [mob_w_func, mob_nw_func, mob_t_func, fw_func, dfw_dsw_func] = ...
    createMultiLayerExtendedFunctions(model_types, param_structs, mu_w, mu_nw)
% CREATEMULTILAYEREXTENDEDFUNCTIONS Create extended function handles for multiple layers
    
    n_layers = length(model_types);
    
    mob_w_func = cell(n_layers,1);
    mob_nw_func = cell(n_layers,1);
    mob_t_func = cell(n_layers,1);
    fw_func = cell(n_layers,1);
    dfw_dsw_func = cell(n_layers,1);


    for layer_idx = 1:n_layers
        % Create mobility functions
        mob_w = @(sw) evaluateLayerMobility(sw, layer_idx,...
            model_types, param_structs, mu_w, mu_nw, 'water');
        mob_w_func{layer_idx} = @(sw) mob_w(sw);
        mob_nw = @(sw) evaluateLayerMobility(sw, layer_idx, ...
            model_types, param_structs, mu_w, mu_nw, 'nonwetting');
         mob_nw_func{layer_idx} =  @(sw) mob_nw(sw);
        mob_t_func{layer_idx} = @(sw) (mob_w(sw) + mob_nw(sw));

        % Create fractional flow function
        fw = @(sw) (mob_w(sw) ./ ...
             (mob_w(sw) + mob_nw(sw) + eps));
        fw_func{layer_idx} = @(sw) fw(sw);

        % Create derivative function using complex step
        dfw_dsw = @(sw) complexStepLayerDerivative(fw, sw);
        dfw_dsw_func{layer_idx} = @(sw) dfw_dsw(sw);
        
    end
    % % Create mobility functions
    % mob_w_func = @(sw, layer_idx) evaluateLayerMobility(sw, layer_idx, ...
    %     model_types, param_structs, mu_w, mu_nw, 'water');
    % mob_nw_func = @(sw, layer_idx) evaluateLayerMobility(sw, layer_idx, ...
    %     model_types, param_structs, mu_w, mu_nw, 'nonwetting');
    % mob_t_func = @(sw, layer_idx) mob_w_func(sw, layer_idx) + mob_nw_func(sw, layer_idx);
    % 
    % % Create fractional flow function
    % fw_func = @(sw, layer_idx) mob_w_func(sw, layer_idx) ./ ...
    %     (mob_w_func(sw, layer_idx) + mob_nw_func(sw, layer_idx) + eps);
    % 
    % % Create derivative function using complex step
    % dfw_dsw_func = @(sw, layer_idx) complexStepLayerDerivative(fw_func, sw, layer_idx);
end

% =========================================================================
% HELPER FUNCTIONS - EVALUATION
% =========================================================================

function krw = evaluateLayerKrw(sw, layer_idx, model_types, param_structs, n_layers)
% EVALUATELAYERKRW Evaluate water relative permeability for specific layer
    
    if nargin < 3 || isempty(layer_idx)
        layer_idx = 1;
    end
    
    if layer_idx < 1 || layer_idx > n_layers
        error('relativePermeability:InvalidLayerIndex', ...
              'Layer index must be between 1 and %d', n_layers);
    end
    
    model_type = standardizeModelType(model_types{layer_idx});
    params = param_structs{layer_idx};
    
    switch model_type
        case 'BROOKS-COREY'
            krw = brooksCoreyKrw(sw, params.lambda, params.Swr, params.Snwr, params.Krw_end);
        case 'VAN-GENUCHTEN'
            krw = vanGenuchtenKrw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krw_end);
        case 'LET'
            krw = letKrw(sw, params.Lw, params.Ew, params.Tw, params.Swr, params.Snwr, params.Krw_end);
        case 'Chierici'
            krw = chiericiKrw(sw, params.A, params.L, params.Swr, params.Snwr, params.Krw_end);
    end
end

function krnw = evaluateLayerKrnw(sw, layer_idx, model_types, param_structs, n_layers)
% EVALUATELAYERKRNW Evaluate non-wetting relative permeability for specific layer
    
    if nargin < 3 || isempty(layer_idx)
        layer_idx = 1;
    end
    
    if layer_idx < 1 || layer_idx > n_layers
        error('relativePermeability:InvalidLayerIndex', ...
              'Layer index must be between 1 and %d', n_layers);
    end
    
    model_type = standardizeModelType(model_types{layer_idx});
    params = param_structs{layer_idx};
    
    switch model_type
        case 'BROOKS-COREY'
            krnw = brooksCoreyKrnw(sw, params.lambda, params.Swr, params.Snwr, params.Krnw_end);
        case 'VAN-GENUCHTEN'
            krnw = vanGenuchtenKrnw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krnw_end);
        case 'LET'
            krnw = letKrnw(sw, params.Lnw, params.Enw, params.Tnw, params.Swr, params.Snwr, params.Krnw_end);
        case 'Chierici'
            krnw = chiericiKrnw(sw, params.B, params.M, params.Swr, params.Snwr, params.Krnw_end);
    end
end

function mobility = evaluateLayerMobility(sw, layer_idx, model_types, param_structs, mu_w, mu_nw, phase)
% EVALUATELAYERMOBILITY Evaluate mobility for specific layer and phase
    
    if nargin < 3 || isempty(layer_idx)
        layer_idx = 1;
    end
    
    model_type = standardizeModelType(model_types{layer_idx});
    params = param_structs{layer_idx};
    
    switch model_type
        case 'BROOKS-COREY'
            if strcmpi(phase, 'water')
                kr = brooksCoreyKrw(sw, params.lambda, params.Swr, params.Snwr, params.Krw_end);
                mobility = kr / mu_w;
            else
                kr = brooksCoreyKrnw(sw, params.lambda, params.Swr, params.Snwr, params.Krnw_end);
                mobility = kr / mu_nw;
            end
            
        case 'VAN-GENUCHTEN'
            if strcmpi(phase, 'water')
                kr = vanGenuchtenKrw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krw_end);
                mobility = kr / mu_w;
            else
                kr = vanGenuchtenKrnw(sw, params.m, params.n, params.Swr, params.Snwr, params.Krnw_end);
                mobility = kr / mu_nw;
            end
            
        case 'LET'
            if strcmpi(phase, 'water')
                kr = letKrw(sw, params.Lw, params.Ew, params.Tw, params.Swr, params.Snwr, params.Krw_end);
                mobility = kr / mu_w;
            else
                kr = letKrnw(sw, params.Lnw, params.Enw, params.Tnw, params.Swr, params.Snwr, params.Krnw_end);
                mobility = kr / mu_nw;
            end
        case 'Chierici'
            if strcmpi(phase, 'water')
                kr = chiericiKrw(sw, params.A, params.L, params.Swr, params.Snwr, params.Krw_end);
                mobility = kr / mu_w;
            else
                kr = chiericiKrnw(sw, params.B, params.M, params.Swr, params.Snwr, params.Krnw_end);
                mobility = kr / mu_nw;
            end
    end
end

% =========================================================================
% NUMERICAL DIFFERENTIATION FUNCTIONS
% =========================================================================

function df = complexStepDerivative(f, x)
% COMPLEXSTEPDERIVATIVE Compute derivative using complex step method
%   More accurate than finite differences, avoids subtractive cancellation
    
    h = 1e-20;  % Complex step size
    
    if isscalar(x)
        df = imag(f(x + h*1i)) / h;
    else
        % Vectorized version
        df = zeros(size(x));
        for i = 1:numel(x)
            xi = x(i);
            df(i) = imag(f(xi + h*1i)) / h;
        end
        df = reshape(df, size(x));
    end
end

function df = complexStepLayerDerivative(f, x)
% COMPLEXSTEPLAYERDERIVATIVE Compute derivative for layer function
    
    h = 1e-20;
    
    if isscalar(x)
        df = imag(f(x + h*1i)) / h;
    else
        df = zeros(size(x));
        for i = 1:numel(x)
            xi = x(i);
            df(i) = imag(f(xi + h*1i)) / h;
        end
        df = reshape(df, size(x));
    end
end

% =========================================================================
% RELATIVE PERMEABILITY MODEL FUNCTIONS (MODIFIED FOR SATURATION RANGE HANDLING)
% =========================================================================

function krw = brooksCoreyKrw(sw, lambda, Swr, Snwr, Krw_end)
% BROOKSCOREYKRW - Brooks-Corey water relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   lambda - Pore size distribution index
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krw_end - Water endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;  % Minimum valid water saturation
sw_max = 1 - Snwr;  % Maximum valid water saturation

% Check for out-of-range saturations and issue warning
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation for sw < sw_min (linear extrapolation to 0)
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end

% Handle extrapolation for sw > sw_max (extrapolate to 1)
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 + (sw(sw > sw_max) - sw_max) / (1 - sw_max);
end

% Ensure Se is within [0, 1] after extrapolation
Se = min(max(Se, 0), 1);

% Compute relative permeability
exponent_w = (2 + 3*lambda) / lambda;
krw_base = Se .^ exponent_w;

% Apply endpoint scaling
krw = Krw_end .* krw_base;

% For extrapolated values above sw_max, ensure krw <= 1
if any(sw(:) > sw_max)
    % Linear extrapolation from endpoint to 1 at sw = 1
    extrap_ratio = min((sw(sw > sw_max) - sw_max) / (1 - sw_max), 1);
    krw(sw > sw_max) = Krw_end + (1 - Krw_end) .* extrap_ratio;
end

% Ensure krw doesn't exceed 1
krw = min(krw, 1);
end

function krnw = brooksCoreyKrnw(sw, lambda, Swr, Snwr, Krnw_end)
% BROOKSCOREYKRNW - Brooks-Corey non-wetting phase relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   lambda - Pore size distribution index
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krnw_end - Non-wetting endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;  % Minimum valid water saturation
sw_max = 1 - Snwr;  % Maximum valid water saturation

% Check for out-of-range saturations and issue warning
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation for sw < sw_min (extrapolate from 1)
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end

% Handle extrapolation for sw > sw_max (extrapolate to 0)
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 ;
end

% Ensure Se is within [0, 1] after extrapolation
Se = min(max(Se, 0), 1);

% Compute relative permeability
exponent_nw = (2 + lambda) / lambda;
krnw_base = (1 - Se).^2 .* (1 - Se.^exponent_nw);

% Apply endpoint scaling
krnw = Krnw_end .* krnw_base;

% For extrapolated values above sw_max, ensure krnw approaches 0
if any(sw(:) > sw_max)
    % Linear extrapolation from endpoint to 0 at sw = 1
    krnw(sw > sw_max) = 0;
end

% For extrapolated values below sw_min, ensure krnw doesn't exceed endpoint
if any(sw(:) < sw_min)
    % Linear extrapolation from 1 at sw = 0 to endpoint at sw_min
    extrap_ratio = min((sw_min - sw(sw < sw_min)) / sw_min, 1);
    krnw(sw < sw_min) = min(Krnw_end + (1 - Krnw_end) .* extrap_ratio, 1);
end

% Ensure krnw doesn't exceed 1
krnw = min(krnw, 1);
end

function krw = vanGenuchtenKrw(sw, m, n, Swr, Snwr, Krw_end)
% VANGENUCHTENKRW - Van Genuchten water relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   m, n - Shape parameters
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krw_end - Water endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 + (sw(sw > sw_max) - sw_max) / (1 - sw_max);
end

Se = min(max(Se, 0), 1);

% Compute relative permeability
term1 = 1 - (1 - Se.^(1/m)).^m;
krw_base = sqrt(Se) .* term1.^2;

% Apply endpoint scaling
krw = Krw_end .* krw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    extrap_ratio = min((sw(sw > sw_max) - sw_max) / (1 - sw_max), 1);
    krw(sw > sw_max) = Krw_end + (1 - Krw_end) .* extrap_ratio;
end

krw = min(krw, 1);
end

function krnw = vanGenuchtenKrnw(sw, m, n, Swr, Snwr, Krnw_end)
% VANGENUCHTENKRNW - Van Genuchten non-wetting phase relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   m, n - Shape parameters
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krnw_end - Non-wetting endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 ;
end

Se = min(max(Se, 0), 1);

% Compute relative permeability
krnw_base = (1 - Se).^(1/3) .* (1 - Se.^(1/m)).^(2*m);

% Apply endpoint scaling
krnw = Krnw_end .* krnw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    krnw(sw > sw_max) = 0;
end
if any(sw(:) < sw_min)
    extrap_ratio = min((sw_min - sw(sw < sw_min)) / sw_min, 1);
    krnw(sw < sw_min) = min(Krnw_end + (1 - Krnw_end) .* extrap_ratio, 1);
end

krnw = min(krnw, 1);
end

function krw = letKrw(sw, Lw, Ew, Tw, Swr, Snwr, Krw_end)
% LETKRW - LET model water relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   Lw, Ew, Tw - LET parameters for water phase
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krw_end - Water endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 + (sw(sw > sw_max) - sw_max) / (1 - sw_max);
end

Se = min(max(Se, 0), 1);

% Compute relative permeability
numerator = Se .^ Lw;
denominator = Se .^ Lw + Ew .* (1 - Se) .^ Tw;
krw_base = numerator ./ (denominator + eps);

% Apply endpoint scaling
krw = Krw_end .* krw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    extrap_ratio = min((sw(sw > sw_max) - sw_max) / (1 - sw_max), 1);
    krw(sw > sw_max) = Krw_end + (1 - Krw_end) .* extrap_ratio;
end

krw = min(krw, 1);
end

function krnw = letKrnw(sw, Lnw, Enw, Tnw, Swr, Snwr, Krnw_end)
% LETKRNW - LET model non-wetting phase relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   Lnw, Enw, Tnw - LET parameters for non-wetting phase
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krnw_end - Non-wetting endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - Swr - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = 1 ;
end

Se = min(max(Se, 0), 1);

% Compute relative permeability
numerator = (1 - Se) .^ Lnw;
denominator = (1 - Se) .^ Lnw + Enw .* Se .^ Tnw;
krnw_base = numerator ./ (denominator + eps);

% Apply endpoint scaling
krnw = Krnw_end .* krnw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    krnw(sw > sw_max) = 0;
end
if any(sw(:) < sw_min)
    extrap_ratio = min((sw_min - sw(sw < sw_min)) / sw_min, 1);
    krnw(sw < sw_min) = min(Krnw_end + (1 - Krnw_end) .* extrap_ratio, 1);
end

krnw = min(krnw, 1);
end

function krw = chiericiKrw(sw, A, L, Swr, Snwr, Krw_end)
% CHIERICIKrw - chierici water relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   A - Normalized saturation multiplier
%   L - Normalized saturation exponent
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krw_end - Water endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - sw - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = inf;
end

Se = max(Se, 0);

% Compute relative permeability
term1 = -A .* (Se.^-L);
krw_base = exp(term1);

% Apply endpoint scaling
krw = Krw_end .* krw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    extrap_ratio = min((sw(sw > sw_max) - sw_max) / (1 - sw_max), 1);
    krw(sw > sw_max) = Krw_end + (1 - Krw_end) .* extrap_ratio;
end

krw = min(krw, 1);
end

function krnw = chiericiKrnw(sw, B, M, Swr, Snwr, Krnw_end)
% CHIERICIKRNW - Chierici non-wetting phase relative permeability
%
% INPUTS:
%   sw - Water saturation (scalar, vector, or matrix)
%   B - Normalized saturation multiplier
%   M - Normalized saturation exponent
%   Swr - Residual water saturation
%   Snwr - Residual non-wetting phase saturation
%   Krnw_end - Non-wetting endpoint relative permeability

% Validate saturation
if min(sw) < 0 || max(sw) > 1
    error('relativePermeability:InputSaturation',...
        'Input saturation must be between 0 and 1');
end

% Calculate valid saturation range boundaries
sw_min = Swr;
sw_max = 1 - Snwr;

% Check for out-of-range saturations
if any(sw(:) < sw_min) || any(sw(:) > sw_max)
    out_of_range = (sw < sw_min) | (sw > sw_max);
    if any(out_of_range(:))
        warning('relativePermeability:SaturationOutOfRange', ...
            'Input saturations outside valid range [%.3f, %.3f]. Extrapolation applied.', ...
            sw_min, sw_max);
    end
end

% Calculate effective saturation for valid range
Se = (min(max(sw, sw_min), sw_max) - Swr) / (1 - sw - Snwr);

% Handle extrapolation
if any(sw(:) < sw_min)
    Se(sw < sw_min) = 0;
end
if any(sw(:) > sw_max)
    Se(sw > sw_max) = inf ;
end

Se = max(Se, 0);

% Compute relative permeability
term1 = -B .* (Se.^M);
krnw_base = exp(term1);

% Apply endpoint scaling
krnw = Krnw_end .* krnw_base;

% Extrapolation handling
if any(sw(:) > sw_max)
    krnw(sw > sw_max) = 0;
end
if any(sw(:) < sw_min)
    extrap_ratio = min((sw_min - sw(sw < sw_min)) / sw_min, 1);
    krnw(sw < sw_min) = min(Krnw_end + (1 - Krnw_end) .* extrap_ratio, 1);
end

krnw = min(krnw, 1);
end
% =========================================================================
% ADDITIONAL HELPER FUNCTIONS
% =========================================================================

function model_type = standardizeModelType(model_type)
% STANDARDIZEMODELTYPE Convert model type to standardized form
    
    model_type = validatestring(upper(model_type), ...
        {'BROOKS-COREY', 'BC', 'VAN-GENUCHTEN', 'VG', 'LET', 'Chierici'});
    
    switch model_type
        case {'BC', 'BROOKS-COREY'}
            model_type = 'BROOKS-COREY';
        case {'VG', 'VAN-GENUCHTEN'}
            model_type = 'VAN-GENUCHTEN';
        case 'LET'
            model_type = 'LET';
        case 'Chierici'
            model_type = 'Chierici';
    end
end

function validateParameters(params, model_type)
% VALIDATEPARAMETERS Validate model parameters
    
    required_fields = {};
    switch model_type
        case 'BROOKS-COREY'
            required_fields = {'lambda', 'Swr', 'Snwr', 'Krw_end', 'Krnw_end'};
        case 'VAN-GENUCHTEN'
            required_fields = {'m', 'n', 'Swr', 'Snwr', 'Krw_end', 'Krnw_end'};
        case 'LET'
            required_fields = {'Lw', 'Ew', 'Tw', 'Lnw', 'Enw', 'Tnw', 'Swr', 'Snwr', 'Krw_end', 'Krnw_end'};
        case 'Chierici'
            required_fields = {'A', 'L', 'B', 'M', 'Swr', 'Snwr', 'Krw_end', 'Krnw_end'};
    end
    
    missing_fields = {};
    for i = 1:length(required_fields)
        if ~isfield(params, required_fields{i})
            missing_fields{end+1} = required_fields{i};
        end
    end
    
    if ~isempty(missing_fields)
        error('relativePermeability:MissingParameters', ...
              'Missing parameters for %s model: %s', ...
              model_type, strjoin(missing_fields, ', '));
    end
    
    % Validate saturation bounds
    if params.Swr < 0 || params.Swr >= 1
        error('relativePermeability:InvalidSwr', ...
              'Swr must be in [0, 1)');
    end
    
    if params.Snwr < 0 || params.Snwr >= 1 || (params.Swr + params.Snwr) >= 1
        error('relativePermeability:InvalidSnwr', ...
              'Snwr must be in [0, 1) and Swr + Snwr < 1');
    end
    
    % Validate endpoint relative permeabilities
    if params.Krw_end < 0 || params.Krw_end > 1
        error('relativePermeability:InvalidKrwEnd', ...
              'Krw_end must be in [0, 1]');
    end
    
    if params.Krnw_end < 0 || params.Krnw_end > 1
        error('relativePermeability:InvalidKrnwEnd', ...
              'Krnw_end must be in [0, 1]');
    end
end

function displayModelInfo(model_types, param_structs, viscosities_provided)
% DISPLAYMODELINFO Display information about configured models
    
    fprintf('\n=== Relative Permeability Models Configuration ===\n');
    fprintf('Number of layers: %d\n', length(model_types));
    
    for i = 1:length(model_types)
        fprintf('\nLayer %d:\n', i);
        fprintf('  Model: %s\n', model_types{i});
        
        params = param_structs{i};
        
        % Display key parameters
        fprintf('  Swr: %.4f\n', params.Swr);
        fprintf('  Snwr: %.4f\n', params.Snwr);
        fprintf('  Valid saturation range: [%.4f, %.4f]\n', ...
            params.Swr, 1 - params.Snwr);
        fprintf('  Krw_end: %.4f\n', params.Krw_end);
        fprintf('  Krnw_end: %.4f\n', params.Krnw_end);
        
        % Display model-specific parameters
        param_names = fieldnames(params);
        for j = 1:length(param_names)
            if ~ismember(param_names{j}, {'Swr', 'Snwr', 'Krw_end', 'Krnw_end'})
                fprintf('  %s: %.4f\n', param_names{j}, params.(param_names{j}));
            end
        end
    end
    
    if viscosities_provided
        fprintf('\nViscosities provided - Extended functions available\n');
    else
        fprintf('\nNote: Provide viscosities to get mobility and fractional flow functions\n');
    end
    fprintf('================================================\n\n');
end