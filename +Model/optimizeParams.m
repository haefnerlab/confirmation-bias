function [optim_params, optim_correct, optim_flag] = ...
    optimizeParams(trials, frames, params, variables)
%OPTIMIZEPARAMS find the best sampling params (in terms of percent correct)
%for the given data-generating params.
%
% params = OPTIMIZEPARAMS(trials, frames, params, variables) searches using
% 'trials' trials per evaluation of correctness. 'variables' is a cell
% array of params fields that may be searched, and defaults to just
% {'p_match'}, but may be any subset of {'p_match', 'var_e', 'gamma',
% 'prior_D'}

if nargin < 4, variables = {'p_match'}; end

    function correct = percent_correct(vars)
        sample_params = params;
        for i=1:length(variables)
            sample_params.(variables{i}) = vars(i);
        end
        % Run sampling model using 'params' to generate data but
        % sample_params for the model itself.
        results = Model.runSamplingModel(...
            Model.genDataWithParams(trials, frames, params), sample_params);
        correct = sum(results.choices == +1) / trials;
    end

% Initialize vars to the generative values.
vars0 = cellfun(@(v) params.(v), variables);

% Run pattern search with lower and upper bounds on variables.
opts = optimoptions(@patternsearch, ...
    'InitialMeshSize', 0.1, ...
    'MaxMeshSize', 0.1, ...
    'CompletePoll', 'on', ...
    'UseParallel', true);
[optim_vars, optim_correct, optim_flag] = patternsearch(...
    @percent_correct, vars0, [], [], [], [], ...
    cellfun(@lower_bound, variables), cellfun(@upper_bound, variables), ...
    opts);

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