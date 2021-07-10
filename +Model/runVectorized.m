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
% - params.lapse = (inference) probability of "lapse" (random choice on a trial regardless of stim)
% - params.temperature = (decision) if > 0, choices are based on raising categorical prob at the end
%   of a trial to 1/temperature. If 0, choices are argmax.
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
% - params.model = 'itb' for "integration to bound" OR 'itb-int' to use the numerically-integrated
%   noise version so that results.prob_choice is a deterministic function of stimuli (useful for
%   model fitting; choice behavior is identical to 'itb' but slower)
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

% Nobody likes dividing by zero; this will scale the resulting LPO by 1/eps, which effectively sets
% all values to +/- infinity, which is the desired behavior; it makes choices equal to the sign of
% LPO.
if params.temperature == 0
    params.temperature = eps;
end

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
    case 'itb-int'
        updateFun = [];
        messageFun = @Model.itbLogOddsMessage;
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
    
    % Map LPO through sigmoid to get choice probability (lapse applied later)
    results.prob_choice = 1./(1+exp(-results.lpo(:,end) / params.temperature));
elseif ~isempty(messageFun)
    % Keep track of a PMF over LPO values per-trial per-frame. Use most extreme values from ideal
    % observer to set range of LPO values.
    extreme_lpo = 3 * max(abs(results.lpo(:)));
    if contains(params.model, 'itb', 'IgnoreCase', true)
        % If this is an ITB model, then there is no sense keeping track of mass beyond the bound. 
        extreme_lpo = min(extreme_lpo, params.bound);
    end
    % Space LPO grid out in increments of about 0.01 (corresponding to a maximum of ~0.25% change in
    % choice probability)
    halfgrid = ceil(extreme_lpo / 0.01);
    % Ensure an even number of edges so that there are an odd number of bins, including exactly 0.
    lpo_grid_edges = linspace(-extreme_lpo, extreme_lpo, 2*halfgrid);
    % Pre- and append bins to accumulate all out-of-bounds probability mass (bound-crossings in ITB
    % models)
    lpo_grid_edges = [-inf lpo_grid_edges +inf];
    % Construct probability mass histogram
    [trials, frames] = size(data);
    % (Maybe) save space by doing LPO updates 'online' over frames
    if params.save_lpo
        results.pr_lpo = zeros(trials, 2*halfgrid+1, frames);
    else
        results.pr_lpo = zeros(trials, 2*halfgrid+1, 1);
    end
    % Initialize to a delta distribution. TODO (?) some width or jitter to starting point if prior_C
    % doesn't fall in the center of a bin?
    log_prior_bin = find(log(prior_C) - log(1 - prior_C) > lpo_grid_edges, 1, 'last');
    results.pr_lpo(:, log_prior_bin, 1) = 1;
    
    % Do propagation of PMF over frames.
    for f=1:frames
        if params.save_lpo
            results.pr_lpo(:, :, f+1) = messageFun(params, data(:, f), lpo_grid_edges, results.pr_lpo(:, :, f));
        else
            results.pr_lpo(:, :, 1) = messageFun(params, data(:, f), lpo_grid_edges, results.pr_lpo(:, :, 1));
        end
    end
    
    % Enforce sanity of PMFs one more time, just in case.
    results.pr_lpo = results.pr_lpo ./ sum(results.pr_lpo, 2);
    
    % Issue warning if unbounded integration but mass is accumulating at the edges.
    if isinf(params.bound) && (any(any(results.pr_lpo(:,1,:) > 1e-3)) || any(any(results.pr_lpo(:,end,:) > 1e-3)))
        warning('Unbounded integrator is hitting edges of discrete space!');
    end
    
    % Map LPO through sigmoid to get choice probability, then sum over PMF (lapse applied later)
    lpo_grid_ctrs = (lpo_grid_edges(1:end-1) + lpo_grid_edges(2:end))/2;
    lpo_grid_ctrs(1) = -params.bound;
    lpo_grid_ctrs(end) = +params.bound;
    
    % Compute and store discrete distribution over the bernoulli 'p' of choice (the pr_lpo field is
    % the distribution of mass, bernoulli_p_grid is its x-axis)
    bernoulli_p_grid = 1./(1+exp(-lpo_grid_ctrs / params.temperature));
    
    % Effective probability of choice is simply the marginal probability across the bernoulli grid
    % (equivalently, think of this in 2 stages: first draw a 'p' from the distribution over 'p's,
    % then draw a choice according to 'p' -- the resulting p(choice) is equivalent to a single draw
    % from the average bernoulli 'p')
    results.prob_choice = sum(results.pr_lpo(:, :, end) .* reshape(bernoulli_p_grid, 1, []), 2);
    
    % To avoid confusion later, remove 'lpo' field which at this point just contains the ideal
    % observer's log odds
    results = rmfield(results, 'lpo');
end

if isfield(params, 'lapse1') && isfield(params, 'lapse2')
    lapse1 = params.lapse1;
    lapse2 = params.lapse2;
else
    lapse1 = params.lapse;
    lapse2 = params.lapse;
end

% Apply lapse: prob choice ranges in [lapse1 lapse2]
results.prob_choice = lapse1 + (1-lapse1-lapse2) * results.prob_choice;
bool_choices = rand(size(results.prob_choice)) < results.prob_choice;
results.choices = sign(double(bool_choices) - .5);
end
