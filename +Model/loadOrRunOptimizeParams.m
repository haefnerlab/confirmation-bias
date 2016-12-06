function [optim_params, optim_correct, optim_flag] = ...
    loadOrRunOptimizeParams(trials, frames, params, variables)

datadir = fullfile('+Model', 'saved results');
if ~exist(datadir, 'dir'), mkdir(datadir); end

prefix = 'optim';

if any(strcmpi('p_match', variables))
    prefix = [prefix '_PM'];
end
if any(strcmpi('var_e', variables))
    prefix = [prefix '_VE'];
end
if any(strcmpi('gamma', variables))
    prefix = [prefix '_G'];
end
if any(strcmpi('prior_D', variables))
    prefix = [prefix '_PD'];
end

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
        Model.optimizeParams(trials, frames, params, variables);
    save(filepath, 'optim_params', 'optim_correct', 'optim_flag');
end
end