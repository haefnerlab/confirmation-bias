function results = runVectorized(params, data)
%RUNVECTORIZED run the model(s) vectorized over trials. Which model to use and its parameters are
%specified in the 'params' struct. Also see Model.newModelParams for more information. The second
%'data' argument, if given, must be a [trials x frames] matrix of signal levels. If 'data' is not
%given, generates new data using Model.genDataWithParams(params), which is repeatable only if
%params.seed is set.
%
%To run the ideal observer:
% - params.

if ~exist('data', 'var')
    % Create a seed so that results.params.seed has a record of the exact data.
    params.seed = randi(1000000000);
    data = Model.genDataWithParams(params);
end

if strcmp(params.model, 'ideal') && params.updates > 1
    params.updates = 1;
end

prior_C = params.prior_C;
updates = params.updates;
noise = params.noise;
gamma = params.gamma;
lapse = params.lapse;

%% Initialize return values

[trials, frames] = size(data);

results = struct(...
    'params', params, ...
    'choices', zeros(trials, 1), ...
    'llo', zeros(trials, frames * updates), ...
    'lpo', zeros(trials, frames * updates + 1));

results.lpo(:, 1) = log(prior_C) - log(1 - prior_C);

switch lower(params.model)
    case 'is'
        logLikeFun = @Model.isLogLikelihood;
    case 'vb'
        logLikeFun = @Model.vbLogLikelihood;
    case 'ideal'
        logLikeOdds = Model.logLikelihoodOdds(params, data);
        results.lpo(:, 2:end) = results.lpo(:, 1) + cumsum(logLikeOdds, 2);
        results.choices = sign(results.lpo(:, end));
        return
    otherwise
        error('Unrecognized model type: %s', params.model);
end

%% Run model (sequentially over frames/updates, but vectorized over trials)

t = 1;
for f=1:frames
    e = data(:, f);
    for u=1:updates
        reults.llo(:, t) = logLikeFun(params, e, results.lpo(:, t));
        % Log-normal random variable with expected value 1
        eta = exp(randn(trials, 1) * noise - noise^2 / 2);
        update_diff = reults.llo(:, t) - gamma * results.lpo(:, t);
        results.lpo(:, t+1) = eta .* (results.lpo(:, t) + update_diff / updates);
        t = t + 1;
    end
end

results.choices = sign(results.lpo(:, end));
if isfield(params, 'lapse')
    lapse_trials = rand(trials, 1) < lapse;
    results.choices(lapse_trials) = rand(sum(lapse_trials), 1) < 0.5;
end
end