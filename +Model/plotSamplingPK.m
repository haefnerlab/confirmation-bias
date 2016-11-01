function plotSamplingPK(trials, frames, p_prior, p_like, gamma, samples, recompute)
%PLOTSAMPLINGPK(trials, frames, p_prior, p_like, gamma, samples) run
%sampling model and plot PK for the given params.

if nargin < 7, recompute = false; end

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 1, trials = 1000; end
if nargin < 2, frames = 12; end
if nargin < 3, p_prior = 0.8; end
if nargin < 4, p_like = 0.8; end

[data, var_e, prefix] = ...
    Model.genDataWithPriorLikelihood(trials, frames, p_prior, p_like);
sampling_params = Model.newSamplingParams(...
    'var_e', var_e, ...
    'p_x', p_prior, ...
    'gamma', gamma, ...
    'samples', samples);

string_id = Model.getModelStringID(prefix, sampling_params);
savename = ['sampling_pk_' string_id '.fig'];
savefile = fullfile(savedir, savename);

if exist(savefile, 'file') && ~recompute
    disp(['Figure exists in ' savename]);
    openfig(savefile);
else
    disp(['Creating new PK and figure in ' savename]);
    [results, data] = Model.loadOrRunSamplingModel(data, prefix, sampling_params);
    
    [weights, ~, errors] = CustomRegression.PsychophysicalKernel(data, results.choices == +1, 1, 0, 5);
    fig = figure();
    errorbar(weights, errors);
    ylim(1.1*[-abs(max(weights)) abs(max(weights))]);
    title(['PK ' strrep(string_id, '_', ' ')]);
    saveas(fig, savefile);
end

end