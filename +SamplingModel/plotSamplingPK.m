function [weights, errors, fig] = plotSamplingPK(params, pk_hprs, ideal_observer, optimize, optim_grid_size)
%PLOTSAMPLINGPK(params) plot PK of sampling model for given params.
%
% [weights, errors, fig] = PLOTSAMPLINGPK(params) returns PK and fig handle
%
% ... = PLOTSAMPLINGPK(params, pk_hprs, ideal, optimize, optim_grid)

savedir = fullfile('+SamplingModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 2, pk_hprs = [1 0 10]; end
if nargin < 3, ideal_observer = false; end
if nargin < 4, optimize = {}; end
if nargin < 5, optim_grid_size = 11; end

optim_prefix = SamplingModel.getOptimPrefix(optimize, optim_grid_size);

results_uid = SamplingModel.getModelStringID(params, ideal_observer);
if isempty(optimize)
    if ~ideal_observer
        results = LoadOrRun(@SamplingModel.runSamplingModelFast, {params}, ...
            fullfile(params.save_dir, results_uid), '-verbose');
    else
        results = LoadOrRun(@SamplingModel.runIdealObserver, {params}, ...
            fullfile(params.save_dir, results_uid), '-verbose');
    end
else
    % Find optimal param settings.
    [optim_params, ~] = LoadOrRun(@SamplingModel.optimizeParams, ...
        {params, optimize, optim_grid_size}, ...
        fullfile(params.save_dir, [optim_prefix '_' results_uid]), ...
        '-verbose');
    % Get model results at the optimal param settings.
    results_uid = SamplingModel.getModelStringID(optim_params);
    results = LoadOrRun(@SamplingModel.runSamplingModelFast, {optim_params}, ...
        fullfile(params.save_dir, results_uid), '-verbose');
end

data = SamplingModel.genDataWithParams(results.params);
[data, choices] = flipTrials(data, results.choices);
regressors = SamplingModel.logLikelihoodOdds(params, data);
[weights, ~, errors] = CustomRegression.PsychophysicalKernel(regressors, choices, ...
    pk_hprs(1), pk_hprs(2), pk_hprs(3));

pk_id = ['PK_' results_uid];
savefile = fullfile(savedir, [pk_id '.fig']);

fig = figure(); hold on;
errorbar(1:params.frames, weights(1:end-1), errors(1:end-1));
errorbar(params.frames+1, weights(end), errors(end), '-r');
xlabel('time');
ylabel('weight');
saveas(fig, savefile);

end

function [data, choices] = flipTrials(data, choices)
flip_indexes = rand(length(choices), 1) < 0.5;
data(flip_indexes, :) = -data(flip_indexes, :);
choices = choices == +1;
choices(flip_indexes) = ~choices(flip_indexes);
end