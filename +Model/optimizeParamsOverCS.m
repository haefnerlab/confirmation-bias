function [optim_params, optim_correct, correct_grid] = optimizeParamsOverCS(category_infos, sensory_infos, params, variables, ngrid)
%OPTIMIZEPARAMSOVERCS like Model.optimizeParams but taking into account smoothness over CS-space

if nargin < 5, ngrid = 21; end

for v=1:length(variables)
    assert(isfield(params, variables{v}), ['params.' variables{v} ' does not exist!']);
end

assert(ngrid > 0, 'BADS over CS space not implemented (yet); use Model.optimizeParams');

% Create grid of 'sensory_infos' and 'category_infos'
[ss, cc] = meshgrid(sensory_infos, category_infos);

% Create grid of param values to search over.
variables_values = cellfun(@(v) linspace(lower_bound(v), upper_bound(v), ngrid), variables, 'UniformOutput', false);
variables_grid = cell(size(variables));
[variables_grid{:}] = ndgrid(variables_values{:});

% Search over all grid points and compute percent correct (really, percent matching the ideal
% observer). Dims are [#sensory #category ngrid ngrid ...] where each optimized variable is one
% 'ngrid' dimension
correct_grid = zeros([size(ss) size(variables_grid{1})]);
parfor idx=1:numel(variables_grid{1}) * numel(ss)
    params_copy = params;
    
    % Get indices into [ss cc] and [variables_grid] separately
    [cs_i, vg_i] = ind2sub([numel(ss) numel(variables_grid{1})], idx);
    
    % Set data-generating parameters.
    params_copy.sensory_info = ss(cs_i);
    params_copy.category_info = cc(cs_i);
    % Set variances for this pair of category- & sensory-info values. (That is, assume that the
    % model 'knows' the environment statistics)
    params_copy.var_s = Model.getEvidenceVariance(ss(cs_i));
    params_copy.p_match = cc(cs_i);
    
    % Look up value for each grid-search variable based on vg_i.
    values = arrayfun(@(v) variables_grid{v}(vg_i), 1:length(variables));
    
    % Get and store results.
    correct_grid(idx) = percent_correct(params_copy, variables, values);
end

% Smooth percent-correct to infer a maximum at each CS-point.
correct_grid_smooth = smoothn(correct_grid);

optim_params = struct();
for v_i=1:length(variables)
    optim_params.(variables{v_i}) = cell(size(ss));
end
optim_correct = zeros(size(ss));
for s_i=1:length(sensory_infos)
    for c_i=1:length(category_infos)
        [optim_correct(s_i, c_i), max_idx] = max(correct_grid_smooth(s_i, c_i, :));
        % Convert from 1d max_idx back to a subscript for each variable.
        max_subs = cell(size(variables));
        [max_subs{:}] = ind2sub(size(variables_grid), max_idx);
        % Look up the value of each variable given its index.
        optim_vars = zeros(size(variables));
        for v=1:length(variables)
            hi = upper_bound(variables{v});
            lo = lower_bound(variables{v});
            % p is the fraction of the range [lo, hi] where the optimal value resides. Must subract
            % 1 in numerator and denominator because matlab is 1-indexed =(
            p = (max_subs{v} - 1) / (ngrid - 1);
            optim_vars(v) = lo + p * (hi - lo);
        end
        
        % Set optimal values in params.
        for j=1:length(variables)
            optim_params(s_i, c_i).(variables{j}) = optim_vars(j);
        end
    end
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
elseif strcmpi(variable, 'var_x')
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
elseif strcmpi(variable, 'var_x')
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
