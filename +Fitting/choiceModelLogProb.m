function [log_post, log_prior, log_lh_per_trial, likelihood_samples] = choiceModelLogProb(params, prior_info, signals, choices, nInner)

if ~iscell(signals)
    signals = {signals};
    choices = {choices};
end

nSets = length(params);
assert(nSets == length(signals));
assert(nSets == length(choices));

likelihood_samples = cell(size(params));
log_lh_per_trial = cell(size(params));

if all(cellfun(@isempty, signals)), return; end

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

%% Evaluate likelihood
log_lh = 0;
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
    
    log_lh_per_trial{iSet} = log(mean(likelihood_samples{iSet}, 2));
    log_lh = log_lh + sum(log_lh_per_trial{iSet});
end

%% Combine prior and likelihood
log_post = log_prior + log_lh;

end