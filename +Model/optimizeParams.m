function [optim_params, optim_correct, optim_flag] = ...
    optimizeParams(trials, frames, params, variables, ngrid)
%OPTIMIZEPARAMS find the best sampling params (in terms of percent correct)
%for the given data-generating params.
%
% params = OPTIMIZEPARAMS(trials, frames, params, variables) searches using
% 'trials' trials per evaluation of correctness. 'variables' is a cell
% array of params fields that may be searched, and defaults to just
% {'p_match'}, but may be any subset of {'p_match', 'var_e', 'gamma',
% 'prior_D'}
%
% OPTIMIZEPARAMS(trials, frames, params, variables, ngrid) performs grid
% search rather than pattern search; each variable is discretized to
% 'ngrid' values between its lower and upper bound.

use_grid = true;

if nargin < 4, variables = {'p_match'}; end
if nargin < 5, use_grid = false; end

    function correct = percent_correct(vars)
        sample_params = params;
        for i=1:length(variables)
            sample_params.(variables{i}) = vars(i);
        end
        % Run sampling model using 'params' to generate data but
        % sample_params for the model itself.
        data = Model.genDataWithParams(trials, frames, params);
        results = Model.runSamplingModel(data, sample_params);
        ideal_results = Model.runIdealObserver(data, params);
        correct = sum(results.choices == ideal_results.choices) / trials;
    end

if use_grid
    % Create grid of param values to search over.
    variables_values = cellfun(@(v) linspace(lower_bound(v), upper_bound(v), ngrid), variables, 'UniformOutput', false);
    variables_grid = cell(size(variables));
    [variables_grid{:}] = ndgrid(variables_values{:});
    % Search over all grid points and compute percent correct (really,
    % percent matching the ideal observer).
    correct_grid = zeros(size(variables_grid{1}));
    for vi=1:numel(variables_grid{1})
        values = arrayfun(@(v) variables_grid{v}(vi), 1:length(variables));
        correct_grid(vi) = percent_correct(values);
    end
    % Smooth and interpolate percent-correct to infer a maximum.
    correct_grid_fine = interpn(smoothn(correct_grid), 3);
    ngrid_fine = length(correct_grid_fine);
    [optim_correct, max_idx] = max(correct_grid_fine(:));
    % Convert from 1d max_idx back to a subscript for each variable.
    max_subs = cell(size(variables));
    [max_subs{:}] = ind2sub(size(correct_grid_fine), max_idx);
    % Linearly look up the value of each variable given its index.
    optim_vars = zeros(size(variables));
    for v=1:length(variables)
        hi = upper_bound(variables{v});
        lo = lower_bound(variables{v});
        % p is the fraction of the range [lo, hi] where the optimal value
        % resides. Must subract 1 in numerator and denominator because
        % matlab is 1-indexed =(
        p = (max_subs{v} - 1) / (ngrid_fine - 1);
        optim_vars(v) = lo + p * (hi - lo);
    end
    optim_flag = ngrid;
else
    % Initialize vars to the generative values.
    vars0 = cellfun(@(v) params.(v), variables);
    
    % Run pattern search with lower and upper bounds on variables.
    opts = optimoptions('patternsearch', ...
        'InitialMeshSize', 0.1, ...
        'MaxMeshSize', 0.1, ...
        'CompletePoll', 'on', ...
        'UseParallel', true);
    [optim_vars, optim_correct, optim_flag] = patternsearch(...
        @percent_correct, vars0, [], [], [], [], ...
        cellfun(@lower_bound, variables), cellfun(@upper_bound, variables), ...
        opts);
end

% Set optimal values in params.
optim_params = params;
for j=1:length(variables)
    optim_params.(variables{j}) = optim_vars(j);
end

end

function lb = lower_bound(variable)
if strcmpi(variable, 'p_match')
    lb = 0.5;
elseif strcmpi(variable, 'var_e')
    lb = 0;
elseif strcmpi(variable, 'gamma')
    lb = 0;
elseif strcmpi(variable, 'prior_D')
    lb = 0;
end
end

function ub = upper_bound(variable)
if strcmpi(variable, 'p_match')
    ub = 1;
elseif strcmpi(variable, 'var_e')
    ub = inf;
elseif strcmpi(variable, 'gamma')
    ub = 1;
elseif strcmpi(variable, 'prior_D')
    ub = 1;
end
end
