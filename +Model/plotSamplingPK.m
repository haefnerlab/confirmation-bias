function [weights, errors, expfit, fig] = plotSamplingPK(trials, frames, params, pk_hprs, ideal_observer, optimize, optim_grid_size)
%PLOTSAMPLINGPK(trials, frames, params, [recompute]) run (or load) sampling
%model and plot PK for the given params.

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 5, ideal_observer = false; end
if nargin < 6, optimize = {}; end
if nargin < 7, optim_grid_size = 11; end

optim_prefix = Model.getOptimPrefix(optimize, optim_grid_size);
[data, data_prefix] = Model.genDataWithParams(trials, frames, params);

string_id = Model.getModelStringID([optim_prefix data_prefix], params);

if isempty(optimize)
    if ~ideal_observer
        [results, data] = Model.loadOrRunSamplingModel(data, data_prefix, params);
    else
        results = Model.runIdealObserver(data, params);
        string_id = ['ideal_' string_id];
    end
else
    optim_params = Model.loadOrRunOptimizeParams(trials, frames, params, optimize, optim_grid_size);
    [results, data] = Model.loadOrRunSamplingModel(data, [optim_prefix, data_prefix], optim_params);
end
[weights, errors, expfit, pk_id] = Model.loadOrRunModelPK(string_id, data, results, pk_hprs);
weights = weights(1:end-1);
errors = errors(1:end-1);

xs = linspace(0, length(weights));
fit = expfit(1) + expfit(2) * exp(-xs / expfit(3));

savefile = fullfile(savedir, [pk_id '.fig']);

fig = figure(); hold on;
plot(xs, fit, '--r', 'LineWidth', 2);
errorbar(weights, errors);
legend('fit', 'weights');
xlabel('time');
ylabel('weight');
ylim(1.1*[-abs(max(weights)+max(errors)) abs(max(weights) + max(errors))]);
title(['PK ' strrep(pk_id, '_', ' ')]);
saveas(fig, savefile);

end