function [sw_front, fw_front, tangent_slope, results] = ...
    findShockFront(Swr, Snwr, fw_func, dfw_dsw_func, SwInitial, varargin)
% FINDSHOCKFRONT Find shock front saturation using bisection method
%
% DESCRIPTION:
%   This function computes the shock front saturation using the Welge tangent
%   method with bisection. It accepts fractional flow and derivative functions
%   directly from the relativePermeability function. Handles both single and
%   multiple layers with comprehensive visualization. The function also
%   handles semi-miscible case described by Noh et al. (2007)
%   <https://onepetro.org/REE/article-abstract/10/04/406/197289/Implications-of-Coupling-Fractional-Flow-and?redirectedFrom=fulltext>. 
%
% ESSENTIAL
%   This funnction uses bisection method as non-linear solver. It's
%   necessary to download and install it from
%   <https://www.mathworks.com/matlabcentral/fileexchange/28150-bisection-method-root-finding>
%
% SYNTAX:
%   For single layer:
%   [sw_front, fw_front, tangent_slope] = findShockFront(Swr, Snwr, fw_func, dfw_dsw_func, SwInitial)
%
%   For multiple layers:
%   [sw_front, fw_front, tangent_slope] = findShockFront(Swr, Snwr, fw_func, dfw_dsw_func, SwInitial, layer_idx)
%
%   With options:
%   [sw_front, fw_front, tangent_slope] = findShockFront(..., 'PropertyName', PropertyValue, ...)
%
% INPUTS:
%   Swr      - Residual water saturation (scalar or vector for multiple layers)
%   Snwr      - Residual non-wettong phase saturation (scalar or vector for multiple layers)
%   fw_func      - Function handle or cell array of function handles for fractional flow
%   dfw_dsw_func - Function handle or cell array for fractional flow derivative
%   SwInitial      - Initial water saturation (scalar or vector for multiple layers)
%
% OPTIONAL INPUTS:
%   layer_idx    - (Optional) Layer index for multi-layer systems (default: 1)
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'LowerBound' - Lower bound for bisection search (scalar or vector for multiple layers)
%   'UpperBound' - Upper bound for bisection search (scalar or vector for multiple layers)
%   'Tolerance'  - Bisection tolerance (default: 1e-8)
%   'NumPoints'  - Number of points for curve discretization that used for plotting (default: 1000)
%   'Process'    -  String varaible determing if process is drainage or
%                   imbibition (Drainage, drainage, Imbibition or imbibition, default: imbibition)
%   'Miscible'    - Activate miscibility and allows starting point of
%                   tanget to be out of saturation and fractional flow
%                   range (true/false, default: false). For miscible case,
%                   initial is visctious point on 45 degree straight
%                   line. This value is estiamted from solubility at
%                   certain conditions
%   'Verbose'    - Display progress (true/false, default: true)
%   'Plot'       - Plot results (true/false or structure with plot options, default: true)
%   'ShowTangent' - Show tangent line on plot (default: true)
%   'ShowFront'   - Show front saturation marker (default: true)
%   'CombinePlots'- Combine multiple layers in single figure with subplots (default: true)
%   'FigureHandle'- Handle to existing figure (default: new figure)
%   'AxesHandles' - Handles to existing axes (default: new axes)
%
% OUTPUTS:
%   sw_front     - Shock front water saturation (scalar for single, vector for multiple)
%   fw_front     - Fractional flow at shock front (scalar for single, vector for multiple)
%   tangent_slope- Slope of tangent line (dfw/dsw at front) (scalar for single, vector for multiple)
%   results      - Structure with detailed results
%
% EXAMPLE:
%   % Single layer
%   params.lambda = 2.0; params.Swr = 0.2; params.Snwr = 0.15;
%   mu_w = 1e-3; mu_nw = 6e-5;
%   [~, ~, ~, ~, ~, fw_func, dfw_dsw_func] = relativePermeability('BC', params, mu_w, mu_nw);
%   [sw_front, fw_front] = findShockFront(params.Swr, params.Snwr, fw_func, dfw_dsw_func, params.Swr, 'Plot', true);
%
%   % Multiple layers with enhanced plot
%   models = {'BC', 'VG'};
%   params{1} = struct('lambda', 2.0, 'Swr', 0.2, 'Snwr', 0.15);
%   params{2} = struct('m', 0.5, 'n', 2.0, 'Swr', 0.25, 'Snwr', 0.10);
%   funcs = relativePermeability(models, params, mu_w, mu_nw);
%   [sw_front, fw_front] = findShockFront([0.2 ; 0.25], [0.15 ; 0.1], funcs.fw_func, funcs.dfw_dsw_func, [0.2 ; 0.25], ...
%       'Plot', struct('ShowTangent', true, 'ShowFront', true, 'CombinePlots', true));
%
% AUTHOR: AbdAllah AbdRabelNaby Youssef ElSaieed
% DATE: 28 December 2025

% =========================================================================
% PARSE INPUTS
% =========================================================================

% Default parameters
defaults.LowerBound = [];
defaults.UpperBound = [];
defaults.Tolerance = 1e-8;
defaults.NumPoints = 1000;
defaults.Process = 'imbibition';
defaults.Miscible = false;
defaults.Verbose = true;
defaults.Plot = true;
defaults.ShowTangent = true;
defaults.ShowFront = true;
defaults.CombinePlots = true;
defaults.FigureHandle = [];
defaults.AxesHandles = [];
defaults.LayerNames = {};

% Check if layer_idx is provided as numeric argument
layer_idx = [];      % initializing layer id with empty vector
if ~isempty(varargin) &&  isscalar(varargin{1}) &&  isPositiveIntegerArray(varargin{1})
    layer_idx = varargin{1};
    varargin = varargin(2:end);
end

% Parse name-value pairs
p = inputParser;
addRequired(p, 'Swr',  @(x) isnumeric(x));
addRequired(p, 'Snwr',  @(x) isnumeric(x));
addRequired(p, 'fw_func');
addRequired(p, 'dfw_dsw_func');
addRequired(p, 'SwInitial',  @(x) isnumeric(x));

fields = fieldnames(defaults);
for i = 1:length(fields)
    addParameter(p, fields{i}, defaults.(fields{i}));
end

parse(p, Swr, Snwr, fw_func, dfw_dsw_func, SwInitial, varargin{:});
params = p.Results;

% Parse plot options if provided as structure
if isstruct(params.Plot)
    plot_opts = params.Plot;
    plot_requested = true;
else
    plot_requested = params.Plot;
    plot_opts = struct();
end



% =========================================================================
% DETERMINE NUMBER OF LAYERS, SETUP AND QC INPUTS
% =========================================================================

% Check if input endpoint saturations are valid
if min(Swr) < 0 || max(Swr) >= 1 || min(Snwr) < 0 || max(Snwr) >= 1 || max(Swr + Snwr) >= 1 || isempty(Swr) || isempty(Snwr)
        error('findShockFront:InvalidResiduals', ...
              ['Input Swr & Snwr must be between zero and 1, ' ...
              '& their sum less than or equal 1']);
end

% Check if miscible option is activated or no
if ~params.Miscible && (max(SwInitial) > 1 || min(SwInitial) < 0)
    error('findShockFront:InvalidInput', ...
        'For immiscible case, initial saturations must be between 0 and 1. Else set "Miscible" to true')
end

% Check if number of points are valid
if ~isPositiveIntegerArray(params.NumPoints) || params.NumPoints < 3
    error('findShockFront:NumPoints', ...
        'Number of points "NumPoints" must be postive integer more than 2');
end

% Check if process is valid
if ~ischar(params.Process)
    error('findShockFront:Process', ...
        '"Process" must be "Drainage" or "Imbibition"');
end
if strcmp(params.Process, 'Drainage') || strcmp(params.Process, 'drainage')
    warning('findShockFront:Process', ...
        'This is draiange process. If you need to switch to imbibition, set "Process" to "imbibition"');
    % Checking input initial saturation for miscible case
    if params.Miscible && min(SwInitial) < 1
        error('findShockFront:InvalidInput', ...
            'For miscible draiange case, initial saturations must be more than 1');
    end
    params.Process = 'drainage';
    UpperBoundCheck = true;         % indicator used to get upper bound later 
    if isempty(params.UpperBound)
        warning('findShockFront:UpperBound',...
            ['Upper bounds are not supplied. Saturation corresponding to maximum fractional flow derivative' ...
            ' will be used instead']);
        UpperBoundCheck = false; 
    end
    LowerBoundCheck = true;         % indicator used to get lower bound later
    if isempty(params.LowerBound)
        warning('findShockFront:LowerBound',...
            'Lower bounds are not supplied. Resiudal water saturation will be used instead');
        LowerBoundCheck = false;    
    end
elseif strcmp(params.Process, 'Imbibition') || strcmp(params.Process, 'imbibition')
    warning('findShockFront:Process', ...
        'This is imbibition process. If you need to switch to draiange, set "Process" to "drainage"');
    % Checking input initial saturation for miscible case
    if params.Miscible && max(SwInitial) > 0
        error('findShockFront:InvalidInput', ...
            'For miscible imbibition case, initial saturations must be less than 0');
    end
    params.Process = 'imbibition';
    LowerBoundCheck = true;         % indicator used to get lower bound later          
    if isempty(params.LowerBound)
        warning('findShockFront:LowerBound',...
            ['Lower bounds are not supplied. Saturation corresponding to maximum fractional flow derivative' ...
            ' will be used instead']);
        LowerBoundCheck = false;    
    end
    UpperBoundCheck = true;         % indicator used to get upper bound later
    if isempty(params.UpperBound)
        warning('findShockFront:UpperBound',...
            ['Upper bounds are not supplied. Saturation corresponding to residual non-wetting phase saturation' ...
            ' will be used instead']);
        UpperBoundCheck = false; 
    end
else 
    error('findShockFront:Process', ...
        '"Process" must be "Drainage" or "Imbibition"');
end

% Check if we have multiple layers
if iscell(fw_func)
    n_layers = length(fw_func);
    
    % Check if the calcualtion of all layers are required
    if isempty(layer_idx)
    % Validate that derivative functions match
        if ~iscell(dfw_dsw_func) || length(dfw_dsw_func) ~= n_layers
            error('findShockFront:InvalidInput', ...
                  'Number of fractional flow functions must match derivative functions');

        % Validate that residual water saturations match and entered as vector
        elseif ~isvector(Swr) || length(Swr) ~= n_layers
            error('findShockFront:InvalidInput', ...
                  'Number of residual water saturations must be a vector match derivative functions');

        % Validate that non-wetting phase residual saturations match and entered as vector
        elseif ~isvector(Snwr) || length(Snwr) ~= n_layers
            error('findShockFront:InvalidInput', ...
                  ['Number of non-wetting phase residual saturations' ...
                  ' must be a vector match derivative functions']);
        
        % Validate that initial saturations match and entered as vector
        elseif ~isvector(SwInitial) || length(SwInitial) ~= n_layers
            error('findShockFront:InvalidInput', ...
                  'Number of initial saturations must be a vector match derivative functions');
       
        % Validate that lower bounds match 
        elseif ~isempty(params.LowerBound)
            if ~isvector(params.LowerBound) || length(params.LowerBound) ~= n_layers
                error('findShockFront:InvalidInput', ...        
                  'Number of lower bounds must be a vector match derivative functions');
            end
    
        % Validate that upper bounds match 
        elseif ~isempty(params.UpperBound)
            if ~isvector(params.UpperBound) || length(params.UpperBound) ~= n_layers
                error('findShockFront:InvalidInput', ...        
                  'Number of upper bounds must match a vector derivative functions');
            end
        end
    
    % Check if calcualtions are required for specific layers
    elseif   ~isPositiveIntegerArray(layer_idx) ||  max(layer_idx) > n_layers
        error('findShockFront:InvalidInput', ...
            'layer idx must be postive integer and less than number of layers')

    else   % only specific layers are required
        n_layers = length(layer_idx);

            % Validate that initial saturations match and entered as vector
            if ~isvector(SwInitial) || length(SwInitial) ~= n_layers
                error('findShockFront:InvalidInput', ...
                      'Number of initial saturations must be a vector match derivative functions');

            % Validate that residual water saturations match and entered as vector
            elseif ~isvector(Swr) || length(Swr) ~= n_layers
                error('findShockFront:InvalidInput', ...
                      'Number of residual water saturations must be a vector match derivative functions');

            % Validate that non-wetting phase residual saturations match and entered as vector
            elseif ~isvector(Snwr) || length(Snwr) ~= n_layers
                   error(['Number of non wetting phase residual saturations' ...
                  ' must be a vector match derivative functions']);
        
            % Validate that lower bounds match 
            elseif ~isempty(params.LowerBound)
                if ~isvector(params.LowerBound) || length(params.LowerBound) ~= n_layers
                    error('findShockFront:InvalidInput', ...        
                      'Number of lower bounds must be a vector match derivative functions');
                end
        
            % Validate that upper bounds match 
            elseif ~isempty(params.UpperBound)
                if ~isvector(params.UpperBound) || length(params.UpperBound) ~= n_layers
                    error('findShockFront:InvalidInput', ...        
                      'Number of upper bounds must match a vector derivative functions');
                end
            end

        fw_func = {fw_func(layer_idx)};
        dfw_dsw_func = {dfw_dsw_func(layer_idx)};

    end

else
    n_layers = 1;
    
    % Convert to cell arrays for uniform handling
    fw_func = {fw_func};
    dfw_dsw_func = {dfw_dsw_func};
end


% Defining tolerance of bisection method 
options.TolX = 0;       % Tolerance in saturation interval
options.TolFun = params.Tolerance;   % Tolerance in target function that captures difference in slope of fractional flow derivative

% Initialize outputs
sw_front = zeros(n_layers, 1);
fw_front = zeros(n_layers, 1);
tangent_slope = zeros(n_layers, 1);

% Prepare results structure
results = struct();
results.n_layers = n_layers;
results.sw_front = zeros(n_layers, 1);
results.fw_front = zeros(n_layers, 1);
results.tangent_slope = zeros(n_layers, 1);
results.converged = false(n_layers, 1);
results.sw_initial = SwInitial;
results.fw_initial = SwInitial;         % This initial value is set a dummy value for immscible case. For miscible case, this is starting point of tangent
results.sw_range = zeros(n_layers, params.NumPoints);
results.fw_range = zeros(n_layers, params.NumPoints);

% =========================================================================
% PROCESS EACH LAYER
% =========================================================================

if params.Verbose
    fprintf('\n=== Shock Front Computation ===\n');
    fprintf('Number of layers: %d\n', n_layers);
    fprintf('Using bisection method with tolerance: %g\n', params.Tolerance);
end


% Looping on layers to get flood flood front saturation
for layer = 1:n_layers
  
    % Get functions for this layer
    current_fw_func = @(sw) fw_func{layer}(sw);
    current_dfw_func = @(sw) dfw_dsw_func{layer}(sw);

    % Get initial saturation for this layer
    sw_initial_layer = SwInitial(layer);
    % Get residual saturations for this layer
    swr_layer = Swr(layer);
    snwr_layer = Snwr(layer);

    if params.Verbose && n_layers > 1
        fprintf('\n--- Processing Layer %d ---\n', layer);
        fprintf('Layer %d - Initial saturation: %.4f - Swr: %.4f - Snwr: %.4f\n', ...
            layer, sw_initial_layer,swr_layer,snwr_layer);
    end


    % Check which process and assign search boundaries
    switch params.Process
        case 'imbibition'
            % Check upper bound
            if ~UpperBoundCheck
                upper_bound_layer = 1 - snwr_layer;
            else
                upper_bound_layer = params.UpperBound(layer);
            end
            % Check lower bound
            if ~LowerBoundCheck
                % Get saturation corresponding to maximum fractional flow
                % derivative and comparing it with initial saturation
                lower_bound_layer = max(fminbnd(@(sw) -1*current_dfw_func(sw),swr_layer,1 - snwr_layer), sw_initial_layer);
            else
                lower_bound_layer = params.LowerBound(layer);
            end
        case 'drainage'
            % Check lower bound
            if ~LowerBoundCheck
                lower_bound_layer = swr_layer;
            else
                lower_bound_layer = params.LowerBound(layer);
            end
            % Check upper bound
            if ~UpperBoundCheck
                % Get saturation corresponding to maximum fractional flow
                % derivative and comparing it with initial saturation
                upper_bound_layer = min(fminbnd(@(sw) -1*current_dfw_func(sw),swr_layer,1 - snwr_layer), sw_initial_layer);
            else
                upper_bound_layer = params.UpperBound(layer);
            end
    end
    


    % =====================================================================
    % BISECTION METHOD FOR TANGENT POINT
    % =====================================================================
    
    if params.Verbose
        fprintf('Layer %d - Bisection bounds: [%.4f, %.4f]\n', layer, lower_bound_layer, upper_bound_layer);
    end
        
    % Define the function to find root of: g(sw) = 0
    g_func = @(sw)  computeTangentCondition(sw, sw_initial_layer, current_fw_func, current_dfw_func,params.Miscible);                                                                            
    
   
    % Use bisection method
    [sw_front_layer, ~, exit_flag] = bisection(...
    g_func, lower_bound_layer, upper_bound_layer, 0 , options);           % if exit flag is 1, 2 or 3, the solution converge
    
    converged = (exit_flag > 0);
    
    if ~converged
        warning('findShockFront:BisectionWarning', ...
                'Layer %d: Bisection did not converge perfectly. Exit flag: %d', layer, exit_flag);
    end
        
    
    % =====================================================================
    % COMPUTE FINAL RESULTS FOR THIS LAYER
    % =====================================================================
    
    % Compute values at shock front
    dfw_front_layer = current_dfw_func(sw_front_layer);
    fw_front_layer = current_fw_func(sw_front_layer);
    sw_range_layer = linspace(swr_layer,1 - snwr_layer,params.NumPoints);
    fw_range_layer = current_fw_func(sw_range_layer);

    % Store results
    sw_front(layer) = sw_front_layer;
    fw_front(layer) = fw_front_layer;
    tangent_slope(layer) = dfw_front_layer;
    
    results.sw_front(layer) = sw_front_layer;
    results.fw_front(layer) = fw_front_layer;
    results.tangent_slope(layer) = dfw_front_layer;
    results.sw_range(layer,:) = sw_range_layer;
    results.fw_range(layer,:) = fw_range_layer;
    if ~params.Miscible
        results.fw_initial(layer) = current_fw_func(sw_initial_layer);
    end


    results.converged(layer) = converged;

    if params.Verbose
        fprintf('Layer %d - Results:\n', layer);
        fprintf('  Shock front saturation: %.6f\n', sw_front_layer);
        fprintf('  Fractional flow at front: %.6f\n', fw_front_layer);
        fprintf('  Tangent slope: %.6f\n', dfw_front_layer);
        fprintf('  Converged: %s\n', string(converged));
    end
end

% =========================================================================
% CREATE COMPREHENSIVE VISUALIZATION
% =========================================================================

if plot_requested
    if params.Verbose
        fprintf('\nGenerating visualization...\n');
    end
    
    % Merge plot options
    final_plot_opts = struct();
    final_plot_opts.ShowTangent = params.ShowTangent;
    final_plot_opts.ShowFront = params.ShowFront;
    final_plot_opts.CombinePlots = params.CombinePlots;
    final_plot_opts.FigureHandle = params.FigureHandle;
    final_plot_opts.AxesHandles = params.AxesHandles;
    final_plot_opts.LayerNames = params.LayerNames;
    
    % Override with user-provided plot options
    if isstruct(params.Plot)
        user_fields = fieldnames(params.Plot);
        for i = 1:length(user_fields)
            final_plot_opts.(user_fields{i}) = params.Plot.(user_fields{i});
        end
    end
    
    % Create visualization
    visualizeShockFronts(results, final_plot_opts);
end

% =========================================================================
% DISPLAY FINAL SUMMARY
% =========================================================================

if params.Verbose
    fprintf('\n=== Final Results ===\n');
    if n_layers == 1
        fprintf('Shock front saturation: %.6f\n', sw_front);
        fprintf('Fractional flow at front: %.6f\n', fw_front);
        fprintf('Tangent slope: %.6f\n', tangent_slope);
    else
        fprintf('Layer\tSw_front\tFw_front\tSlope\t\tConverged\n');
        fprintf('-----\t--------\t--------\t-----\t\t--------\n');
        for layer = 1:n_layers
            fprintf('%d\t%.6f\t%.6f\t%.6f\t%s\t\t\n', ...
                layer, sw_front(layer), fw_front(layer), ...
                tangent_slope(layer), string(results.converged(layer)));
        end
    end
    fprintf('========================\n\n');
end

end % Main function


% =========================================================================
% HELPER FUNCTIONS - BISECTION AND SEARCH
% =========================================================================

function g = computeTangentCondition(sw, sw_initial, fw, dfw_dsw, Miscible)
% COMPUTETANGENTCONDITION Compute tangent condition function

    fw_fun = fw(sw);
    dfw_dsw_fun = dfw_dsw(sw);
    if ~Miscible
        fw_initial = fw(sw_initial);
    else
        fw_initial = sw_initial;   % point on 45 degree line with coordiante (sw_initial , sw_initial); See Noh et al. (2007)
    end
    

    chord_slope = (fw_initial - fw_fun) ./ (sw_initial - sw);
    g = chord_slope - dfw_dsw_fun;
end


% =========================================================================
% VISUALIZATION FUNCTIONS
% =========================================================================

function visualizeShockFronts(results, plot_opts)
% VISUALIZESHOCKFRONTS Create comprehensive visualization
    
    n_layers = results.n_layers;
    
    % Determine figure layout
    if plot_opts.CombinePlots && n_layers > 1
        % Combined figure with subplots
        createCombinedFigure(results, plot_opts);
    else
        % Separate figures for each layer
        for layer = 1:n_layers
            createSingleLayerFigure(results, layer, plot_opts);
        end
    end
end

function createCombinedFigure(results, plot_opts)
% CREATECOMBINEDFIGURE Create figure with subplots for multiple layers
    
    n_layers = results.n_layers;
    
    % Determine subplot layout
    n_cols = min(3, n_layers);
    n_rows = ceil(n_layers / n_cols);
    
    % Create or use existing figure
    if isempty(plot_opts.FigureHandle) || ~isvalid(plot_opts.FigureHandle)
        fig = figure('Position', [100, 100, 1600, 900], 'Color', 'white');
    else
        fig = plot_opts.FigureHandle;
        figure(fig);
        clf(fig);
    end
    
    % Create title
    sgtitle('Shock Front Analysis - Multiple Layers', 'FontSize', 16, 'FontWeight', 'bold');
    
    for layer = 1:n_layers
        % Create or use existing axes
        if isempty(plot_opts.AxesHandles) || length(plot_opts.AxesHandles) < layer
            ax = subplot(n_rows, n_cols, layer);
        else
            ax = plot_opts.AxesHandles(layer);
            axes(ax);
        end
        
        % Plot layer results
        plotLayerResults(ax, results,  layer, plot_opts);
        
        % Add layer title
        if ~isempty(plot_opts.LayerNames) && length(plot_opts.LayerNames) >= layer
            layer_name = plot_opts.LayerNames{layer};
        else
            layer_name = sprintf('Layer %d', layer);
        end
        
        title(ax, sprintf('%s\nsw_{front} = %.4f, fw_{front} = %.4f', ...
            layer_name, results.sw_front(layer), results.fw_front(layer)), ...
            'FontSize', 12);
    end
    
    % Add overall legend to first subplot
    axes(subplot(n_rows, n_cols, 1));
    legend('show', 'Location', 'best');
    
    % Add summary text box
    annotation_str = createSummaryAnnotation(results);
    annotation(fig, 'textbox', [0.02, 0.02, 0.3, 0.1], ...
        'String', annotation_str, 'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 9);
end

function createSingleLayerFigure(results, layer, plot_opts)
% CREATESINGLELAYERFIGURE Create figure for single layer
    
    % Create or use existing figure
    if isempty(plot_opts.FigureHandle) || ~isvalid(plot_opts.FigureHandle)
        fig = figure('Position', [100, 100, 1200, 800], 'Color', 'white');
    else
        fig = plot_opts.FigureHandle;
        figure(fig);
    end
    
    % Create or use existing axes
    if isempty(plot_opts.AxesHandles) || length(plot_opts.AxesHandles) < layer
        ax = axes('Parent', fig);
    else
        ax = plot_opts.AxesHandles(layer);
        axes(ax);
    end
    
    % Plot layer results
    plotLayerResults(ax, results,  layer, plot_opts);
    
    % Add title
    if ~isempty(plot_opts.LayerNames) && length(plot_opts.LayerNames) >= layer
        layer_name = plot_opts.LayerNames{layer};
    else
        layer_name = sprintf('Layer %d', layer);
    end
    
    title(ax, sprintf('%s', layer_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Add results annotation
    annotation_str = sprintf(['Results:\n' ...
                             'sw_{front} = %.4f\n' ...
                             'fw_{front} = %.4f\n' ...
                             'dfw/dsw = %.4f\n' ...
                             'Converged: %s'], ...
                            results.sw_front(layer), results.fw_front(layer), ...
                            results.tangent_slope(layer),  string(results.converged(layer)));
    
    annotation(fig, 'textbox', [0.02, 0.02, 0.25, 0.15], ...
        'String', annotation_str, 'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
end

function plotLayerResults(ax, results,  layer, plot_opts)
% PLOTLAYERRESULTS Plot results for a specific layer
    
    sw_range = results.sw_range(layer,:);
    fw_range = results.fw_range(layer,:);
    % dfw_dsw_range = results.dfw_dsw_range{layer};
    sw_front = results.sw_front(layer);
    fw_front = results.fw_front(layer);
    sw_initial = results.sw_initial(layer);
    fw_initial = results.fw_initial(layer);
    
    % Plot fractional flow curve
    plot(ax, sw_range, fw_range, 'm-', 'LineWidth', 3, 'DisplayName', 'f_w');
    hold(ax, 'on');
    
    % % Plot derivative curve
    % yyaxis(ax, 'right');
    % plot(ax, sw_range, dfw_dsw_range, 'm--', 'LineWidth', 1.5, ...
    %     'DisplayName', 'df_w/ds_w');
    % ylabel(ax, 'df_w/ds_w');
    % 
    % % Set right y-axis color
    % ax.YAxis(2).Color = [0.5, 0, 0.5];
    
    % Return to left y-axis
    % yyaxis(ax, 'left');
    ylabel(ax, 'f_w');
    
    % Plot initial condition
    plot(ax, sw_initial, fw_initial, 'bs', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Initial (s_{w0})');
    
    % Plot shock front
    if plot_opts.ShowFront
        plot(ax, sw_front, fw_front, 'ro', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'r', 'LineWidth', 1.5, ...
            'DisplayName', 'Shock Front');
    end
    
    % Plot tangent line
    if plot_opts.ShowTangent
        slope = results.tangent_slope(layer);
        % x_tangent = linspace(sw_front, sw_initial, 50);
        x_tangent = linspace(0, 1, 50);
        y_tangent = fw_front - slope * (sw_front - x_tangent);
        
        plot(ax, x_tangent, y_tangent, 'k--', 'LineWidth', 2.5, ...
            'DisplayName', 'Tangent Line');
    end
    


    % Formatting
    xlabel(ax, 's_w');
    grid(ax, 'on');
    box(ax, 'on');
    
    % Set limits
    xlim(ax, [0 , 1]);
    ylim_left = [0, 1];
    ylim(ax, ylim_left);
    
    % % Add second y-axis limits
    % yyaxis(ax, 'right');
    % y_deriv = dfw_dsw_range;
    % ylim(ax, [min(y_deriv)*0.9, max(y_deriv)*1.1]);
    
    % Return to left y-axis for legend
    % yyaxis(ax, 'left');
    legend(ax, 'Location', 'best', 'FontSize', 15);
    set(gca,'FontSize',20,'FontWeight','bold','linewidth',3); set(gcf,'color','w');
    hold(ax, 'off');
end

function annotation_str = createSummaryAnnotation(results)
% CREATESUMMARYANNOTATION Create summary annotation text
    
    n_layers = results.n_layers;
    
    if n_layers == 1
        annotation_str = sprintf(['Single Layer Results:\n' ...
                                 'sw_{front} = %.4f\n' ...
                                 'fw_{front} = %.4f\n' ...
                                 'Converged: %s'], ...
                                results.sw_front(1), results.fw_front(1), ...
                                string(results.converged(1)));
    else
        annotation_str = 'Multi-Layer Summary:\n';
        for layer = 1:min(5, n_layers)  % Show first 5 layers in summary
            annotation_str = sprintf('%sL%d: sw=%.4f, fw=%.4f\n', ...
                annotation_str, layer, results.sw_front(layer), results.fw_front(layer));
        end
        
        if n_layers > 5
            annotation_str = sprintf('%s... and %d more layers', annotation_str, n_layers - 5);
        end
        
        % Add convergence summary
        n_converged = sum(results.converged);
        annotation_str = sprintf('%s\nConverged: %d/%d layers', ...
            annotation_str, n_converged, n_layers);
    end
end