function [optim_params, optim_correct, optim_flag, prefix] = ...
    loadOrRunOptimizeParams(trials, frames, params, variables, grid_search_size)

datadir = fullfile('+Model', 'saved results');
if ~exist(datadir, 'dir'), mkdir(datadir); end

if nargin < 5, grid_search_size = []; end

prefix = Model.getOptimPrefix(variables, grid_search_size);
string_id = Model.getModelStringID(prefix, params);

filename = [string_id '.mat'];
filepath = fullfile(datadir, filename);

if exist(filepath, 'file')
    disp(['Loading precomputed results for ' filename]);
    contents = load(filepath);
    optim_params = contents.optim_params;
    optim_correct = contents.optim_correct;
    optim_flag = contents.optim_flag;
else
    disp(['Computing new results for ' filename]);
    [optim_params, optim_correct, optim_flag] = ...
        Model.optimizeParams(trials, frames, params, variables, grid_search_size);
    save(filepath, 'optim_params', 'optim_correct', 'optim_flag');
end
end