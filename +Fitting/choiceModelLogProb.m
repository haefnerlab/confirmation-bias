function [log_post, log_like, est_variance, likelihood_samples] = choiceModelLogProb(params, prior_info, signals, choices, nInner)

if ~iscell(signals)
    signals = {signals};
    choices = {choices};
end

nSets = length(params);
assert(nSets == length(signals));
assert(nSets == length(choices));

likelihood_samples = cell(size(params));
log_like_per_trial = cell(size(params));

if all(cellfun(@isempty, signals)), return; end

if Model.isStochastic(params) && nInner == 1
    warning('It''s strongly suggested that # inner loop simulations > 1 when using a stochastic model!');
end

% Sanity check that all choices are in [-1 +1]
for iSet=1:nSets
    if ~all(ismember(choices{iSet}, [-1 +1]))
        error('Choices in wrong space! Should be in [-1 +1] but seeing %s.', mat2str(unique(choices{iSet})));
    end
end

%% Evaluate prior
log_prior = 0;
prior_fields = fieldnames(prior_info);
for iF=1:length(prior_fields)
    field = prior_fields{iF};
    val = Fitting.getParamsFields(params(1), field);
    if isfield(prior_info.(field), 'logpriorpdf')
        log_prior = log_prior + prior_info.(field).logpriorpdf(val);
    else
        log_prior = log_prior + log(prior_info.(field).priorpdf(val));
    end
end

%% Evaluate likelihood per trial
for iSet=1:nSets
    thisParams = params(iSet);
    nTrials = length(choices{iSet});
    
    % Repeat signals to vectorize 'nInner' as if it were more trials
    signals{iSet} = repmat(signals{iSet}, nInner, 1);
    
    sim_results = Model.runVectorized(thisParams, signals{iSet});
    if isfield(thisParams, 'lapse_1')
        lapse1 = thisParams.lapse_1;
        lapse2 = thisParams.lapse_2;
    else
        lapse1 = thisParams.lapse;
        lapse2 = thisParams.lapse;
    end

    prob_choice = lapse1 + (1-lapse1-lapse2) ./ (1 + exp(-sim_results.lpo(:, end) / thisParams.temperature));
    
    % Reshape prob_choice from model to [nTrials x nInner]
    prob_choice = reshape(prob_choice, nTrials, nInner);
    
    % likelihood_samples is [nTrials x nInner] and contains the full choice likelihood for each run of
    % the model on each choice.
    chosePos = choices{iSet} == +1;
    choseNeg = ~chosePos;
    likelihood_samples{iSet}(chosePos, :) = prob_choice(chosePos, :);
    likelihood_samples{iSet}(choseNeg, :) = 1-prob_choice(choseNeg, :);
end

% Correct for log transformation. Above, each 'likelihood_samples' is an unbiased estimate of the
% likelihood. However, we want an estimate of the log likelihood. With some assumptions, the
% following corrects for the log transform (see https://stats.stackexchange.com/q/57766). If
% nInner=1, then nothing happens since var_like_per_trial will be 0; this should only be the case if
% using a deterministic model (otherwise, a warning will be issued above).
all_likelihood_samples = vertcat(likelihood_samples{:});
mean_like_per_trial = mean(all_likelihood_samples, 2);
var_like_per_trial = var(all_likelihood_samples, [], 2) / nInner;
est_log_like_per_trial = log(mean_like_per_trial); %- var_like_per_trial./(2*mean_like_per_trial.^2);
est_var_log_like_per_trial = var_like_per_trial./mean_like_per_trial.^2;

log_like = sum(est_log_like_per_trial);
est_variance = sum(est_var_log_like_per_trial);

%% Combine prior and likelihood
log_post = log_prior + log_like;

end