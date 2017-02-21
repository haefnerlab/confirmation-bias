function [optim_params, optim_correct, optim_flag, prefix] = ...
    loadOrRunOptimizeParams(trials, frames, params, variables, grid_search_size)

datadir = fullfile('+Model', 'saved results');
if ~exist(datadir, 'dir'), mkdir(datadir); end

prefix = 'optim';

if exist('grid_search_size', 'var')
    prefix = [prefix '_grid' num2str(grid_search_size)];
end

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
    if exist('grid_search_size', 'var')
        [optim_params, optim_correct, optim_flag] = ...
            Model.optimizeParams(trials, frames, params, variables, grid_search_size);
    else
        [optim_params, optim_correct, optim_flag] = ...
            Model.optimizeParams(trials, frames, params, variables);
    end
    save(filepath, 'optim_params', 'optim_correct', 'optim_flag');
end
end