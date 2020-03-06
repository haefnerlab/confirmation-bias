function [weights, errors, paramfit, fig] = plotPK(params, pk_hprs, optimize, optim_grid_size, do_plot)
%PLOTPK(params) plot PK of sampling model for given params.
%
% [weights, errors, fig] = PLOTPK(params) returns PK and fig handle
%
% ... = PLOTPK(params, pk_hprs, ideal, optimize, optim_grid)

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 2, pk_hprs = [0 0 0]; end
if nargin < 3, optimize = {}; end
if nargin < 4, optim_grid_size = 11; end
if nargin < 5, do_plot = true; end

optim_prefix = Model.getOptimPrefix(optimize, optim_grid_size);

results_uid = Model.getModelStringID(params);
if isempty(optimize)
    results = LoadOrRun(@Model.runVectorized, {params}, ...
        fullfile(params.save_dir, results_uid), '-verbose');
else
    % Find optimal param settings.
    [optim_params, ~] = LoadOrRun(@Model.optimizeParams, ...
        {params, optimize, optim_grid_size}, ...
        fullfile(params.save_dir, [optim_prefix '_' results_uid]), '-verbose');
    % Get model results at the optimal param settings.
    results_uid = Model.getModelStringID(optim_params);
    results = LoadOrRun(@Model.runVectorized, {optim_params}, ...
        fullfile(params.save_dir, results_uid), '-verbose');
end

% Regenerate data from seed
data = Model.genDataWithParams(results.params);
regressors = Model.logLikelihoodOdds(params, data);

% Do regression
if ischar(pk_hprs) && startsWith(pk_hprs, 'exp')
    % TODO - get errors in weight space?
    abb = CustomRegression.ExponentialPK(regressors, results.choices==+1);
    weights = abb(1) * exp(abb(2) * (0:params.frames-1));
    weights(end+1) = abb(3);
    errors = nan(size(weights));
    paramfit = abb;
elseif ischar(pk_hprs) && startsWith(pk_hprs, 'lin')
    % TODO - get errors in weight space?
    sob = CustomRegression.LinearPK(regressors, results.choices==+1);
    weights = sob(2) + sob(1)*(0:params.frames-1);
    weights(end+1) = sob(3);
    errors = nan(size(weights));
    paramfit = sob;
else
    [weights, ~, errors] = CustomRegression.PsychophysicalKernel(regressors, results.choices==+1, ...
        pk_hprs(1), pk_hprs(2), pk_hprs(3));
    paramfit = [];
end

if do_plot || nargout >= 4
    pk_id = ['PK_' results_uid];
    savefile = fullfile(savedir, [pk_id '.fig']);
    fig = figure(); hold on;
    errorbar(1:params.frames, weights(1:end-1), errors(1:end-1));
    errorbar(params.frames+1, weights(end), errors(end), '-r');
    xlabel('time');
    ylabel('weight');
    saveas(fig, savefile);
end

end
