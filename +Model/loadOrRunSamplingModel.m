function [results, data] = loadOrRunSamplingModel(data, prefix, sampling_params, recompute)
%LOADORRUNSAMPLINGMODEL wrapper around runSamplingModel that will load pre-
%existing results if they are available or compute them from the given data
%
% [results, data] = LOADORRUNSAMPLINGMODEL(data, prefix, sampling_params)
% First looks for existing data from the unique combination of the given
% 'prefix' string and 'sampling_params'. Computes new results if there are
% no existing results. If there are precomputed results, the exact data
% used to compute them are returned as well (the 'data' return value);
% 'prefix' should be a string that describes the data.
%
% [results, data] = LOADORRUNSAMPLINGMODEL(..., recompute) if set to true,
% forces a new run of the sampling model.
if nargin < 4, recompute = false; end

datadir = fullfile('+Model', 'saved results');
if ~exist(datadir, 'dir'), mkdir(datadir); end

savename = ['sampling_run_' Model.getModelStringID(prefix, sampling_params) '.mat'];
savefile = fullfile(datadir, savename);

if exist(savefile, 'file') && ~recompute
    disp(['Loading precomputed results from ' savename]);
    contents = load(savefile);
    data = contents.data;
    results = contents.results;
else
    disp(['Computing new results for ' savename]);
    results = Model.runSamplingModel(data, sampling_params);
    save(savefile, 'data', 'results');
end

end