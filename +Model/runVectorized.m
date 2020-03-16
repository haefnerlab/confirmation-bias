function results = runVectorized(params, data)
%RUNVECTORIZED run the model(s) vectorized over trials. Which model to use and its parameters are
%specified in the 'params' struct. Also see Model.newModelParams for more information. The second
%'data' argument, if given, must be a [trials x frames] matrix of signal levels. If 'data' is not
%given, generates new data using Model.genDataWithParams(params), which is repeatable only if
%params.seed is set.
%
%To run the ideal observer:
% - params.model = 'ideal'
% - params.category_info = (data gen) category info
% - params.sensory_info = (data gen) sensory info
% - params.p_match = (inference) category info
% - params.var_s + params.var_x = (inference) sensory info (see @Model.getEvidenceVariance)
% - params.prior_C = (inference) prior probability that C=+1
%
%To run the sampling model:
% - params.model = 'is' for "importance sampling"
% - params.category_info = (data gen) category info
% - params.sensory_info = (data gen) sensory info
% - params.p_match = (inference) category info
% - params.var_x = (inference) variance of each mode of x|C Mixture of Gaussians
% - params.var_s = (inference) variance of s|x
% - params.updates = (inference) number of updates per frame ("n_U" in the paper)
% - params.samples = (inference) number of parallel samples in a single update ("S" in the paper)
% - params.gamma = (inference) linear approximate bias correction term (i.e. leak term)
% - params.noise = (inference) log-normal noise per update is exp(randn*noise - noise^2/2)
% - params.lapse = (inference) probability of "lapse" (random choice on a trial regardless of stim)
% - params.temperature = (decision) if > 0, choices are based on raising categorical prob at the end
%   of a trial to 1/temperature. If 0, choices are argmax.
%
%To run the variational model:
% - params.model = 'vb' for "variational bayes" using the q(x,z)q(C) factorization, or 'vb-cxz' to
%   use the q(x)q(z)q(C) factorization
% - params.category_info = (data gen) category info
% - params.sensory_info = (data gen) sensory info
% - params.p_match = (inference) category info
% - params.var_x = (inference) variance of each mode of x|C Mixture of Gaussians
% - params.var_s = (inference) variance of s|x
% - params.updates = (inference) number of coordinate ascent steps per frame ("n_U" in the paper)
% - params.gamma = (inference) linear approximate bias correction term (i.e. leak term)
% - params.noise = (inference) log-normal noise per update is exp(randn*noise - noise^2/2)
% - params.lapse = (inference) probability of "lapse" (random choice on a trial regardless of stim)
% - params.step_size = (inference) fraction of full update to log odds C to apply in each iteration.
% - params.temperature = (decision) if > 0, choices are based on raising categorical prob at the end
%   of a trial to 1/temperature. If 0, choices are argmax.
%
%To run the ITB model:
% - params.model = 'itb' for "integration to bound"
% - params.category_info = (data gen) category info
% - params.sensory_info = (data gen) sensory info
% - params.p_match = (inference) category info
% - params.var_x = (inference) variance of each mode of x|C Mixture of Gaussians
% - params.var_s = (inference) variance of s|x
% - params.gamma = (inference) linear approximate bias correction term (i.e. leak term)
% - params.noise = (inference) log-normal noise per update is exp(randn*noise - noise^2/2)
% - params.lapse = (inference) probability of "lapse" (random choice on a trial regardless of stim)
% - params.bound = (inference) sticky bound on log posterior odds. Once LPO crosses +/- bound, it
%   stays there to the end of the trial regardless of subsequent frames.
% - params.temperature = (decision) if > 0, choices are based on raising categorical prob at the end
%   of a trial to 1/temperature. If 0, choices are argmax.

if ~exist('data', 'var')
    if ~isfield(params, 'seed') || isempty(params.seed)
        % Create a seed so that results.params.seed has a record of the exact data, unless
        % params.seed is already set
        params.seed = randi(1000000000);
    end
    data = Model.genDataWithParams(params);
else
    % Data matrix is given, so seed is unknown (except perhaps by the calling function, but we don't
    % want results.params.seed to be misleading)
    params.seed = [];
end

prior_C = params.prior_C;
lapse = params.lapse;

%% Initialize return values

[trials, frames] = size(data);

results = struct(...
    'params', params, ...
    'choices', zeros(trials, 1), ...
    'ideal_choices', zeros(trials, 1), ...
    'lpo', zeros(trials, frames + 1));

results.lpo(:, 1) = log(prior_C) - log(1 - prior_C);

%% Compute ideal observer results (fast)

logLikeOdds = Model.logLikelihoodOdds(params, data);
results.lpo(:, 2:end) = results.lpo(:, 1) + cumsum(logLikeOdds, 2);
results.ideal_choices = sign(results.lpo(:, end));

%% Handle different models

switch lower(params.model)
    case 'is'
        updateFun = @Model.isLogOddsUpdate;
    case 'vb'
        updateFun = @Model.vbLogOddsUpdate;
    case 'vb-czx'
        updateFun = @Model.vbLogOddsUpdateCZX;
    case 'itb'
        updateFun = @Model.itbLogOddsUpdate;
    case 'ideal'
        % Nothing to do - use 'lpo' from above
        updateFun = [];
    otherwise
        error('Unrecognized model type: %s', params.model);
end

%% Run model (sequentially over frames/updates, but vectorized over trials)

if ~isempty(updateFun)
    for f=1:frames
        results.lpo(:, f+1) = updateFun(params, data(:, f), results.lpo(:, f));
    end
end

if params.temperature == 0
    results.choices = sign(results.lpo(:, end));
else
    adjusted_lpo = results.lpo(:, end) / params.temperature;
    prob_choice = 1./(1+exp(-adjusted_lpo));
    bool_choices = rand(size(prob_choice)) < prob_choice;
    results.choices = sign(double(bool_choices) - .5);
end
fifty_fifty = results.lpo(:,end) == 0;
random_trials = fifty_fifty | rand(trials, 1) < lapse;
results.choices(random_trials) = sign(rand(sum(random_trials), 1) - 0.5);
end
