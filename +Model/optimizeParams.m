function [optim_params, optim_correct, correct_grid] = optimizeParams(params, variables, ngrid)
%MODEL.OPTIMIZEPARAMS find the best sampling params (in terms of percent correct) for the
%given data-generating params.
%
% params = MODEL.OPTIMIZEPARAMS(params, variables, ngrid) use 'variables' as a cell array of params
% fields that may be searched,; it may be any subset of {'p_match', 'var_s', 'gamma', 'prior_C'}.
% Each variable is discretized into a grid of 'ngrid' equally spaced values. If ngrid is 0, attempts
% to use the BADS algorithm.

for v=1:length(variables)
    assert(isfield(params, variables{v}), ['params.' variables{v} ' does not exist!']);
end

if ngrid > 0
    % Create grid of param values to search over.
    variables_values = cellfun(@(v) linspace(lower_bound(v), upper_bound(v), ngrid), variables, 'UniformOutput', false);
    variables_grid = cell(size(variables));
    [variables_grid{:}] = ndgrid(variables_values{:});
    % Search over all grid points and compute percent correct (really, percent
    % matching the ideal observer).
    correct_grid = zeros(size(variables_grid{1}));
    parfor vi=1:numel(variables_grid{1})
        values = arrayfun(@(v) variables_grid{v}(vi), 1:length(variables));
        correct_grid(vi) = percent_correct(params, variables, values);
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
        % p is the fraction of the range [lo, hi] where the optimal value resides. Must subract 1 in
        % numerator and denominator because matlab is 1-indexed =(
        p = (max_subs{v} - 1) / (ngrid_fine - 1);
        optim_vars(v) = lo + p * (hi - lo);
    end
else
    % Initialize vars to the generative values.
    vars0 = cellfun(@(v) params.(v), variables);
    
    pcwrapper = @(vals) percent_correct(params, variables, vals);
    
    % Run BADS with lower and upper bounds on variables.
    opts = bads('defaults');
    opts.UncertaintyHandling = 1;
    % bernoulli trials with p=.5 have variance .25. Divide variance by # indpendent trials, then
    % take sqrt for noise 'sigma'
    opts.NoiseSize = sqrt(.25 / params.trials);
    [optim_vars, optim_correct] = bads(pcwrapper, vars0, ...
        cellfun(@lower_bound, variables), cellfun(@upper_bound, variables), [], [], [], opts);
end

% Set optimal values in params.
optim_params = params;
for j=1:length(variables)
    optim_params.(variables{j}) = optim_vars(j);
end

end

function correct = percent_correct(params, variables, values)
% Run sampling model on current setting of variables.
for i=1:length(variables)
    params.(variables{i}) = values(i);
end

% TODO - smarter resetting of seed
params.seed = randi(1000000000);

results_uid = Model.getModelStringID(params);
results = LoadOrRun(@Model.runVectorized, {params}, fullfile(params.save_dir, results_uid));

% Rather than compute % correct (which is noisy), compute % agreement with ideal observer (which is
% more rubust at low information levels).
correct = mean(results.choices == results.ideal_choices);
end

function lb = lower_bound(variable)
if strcmpi(variable, 'p_match')
    lb = 0.5;
elseif strcmpi(variable, 'var_s')
    lb = 0;
elseif strcmpi(variable, 'gamma')
    lb = 0;
elseif strcmpi(variable, 'prior_C')
    lb = 0;
elseif strcmpi(variable, 'noise')
    lb = 0;
elseif strcmpi(variable, 'step_size')
    lb = 0;
end
end

function ub = upper_bound(variable)
if strcmpi(variable, 'p_match')
    ub = 1;
elseif strcmpi(variable, 'var_s')
    ub = 10;
elseif strcmpi(variable, 'gamma')
    ub = 1;
elseif strcmpi(variable, 'prior_C')
    ub = 1;
elseif strcmpi(variable, 'noise')
    ub = 10;
elseif strcmpi(variable, 'step_size')
    ub = 1;
end
end