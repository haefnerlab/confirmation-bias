function [log_likelihood, log_likelihood_per_trial, likelihood_samples] = choiceModelLogLikelihood(params, signals, choices, nInner)

if ~iscell(signals)
    signals = {signals};
    choices = {choices};
end

assert(length(params) == length(signals));
assert(length(params) == length(choices));

likelihood_samples = cell(size(params));
log_likelihood_per_trial = cell(size(params));
log_likelihood = 0;

if all(cellfun(@isempty, signals)), return; end

for iSet=1:length(params)
    thisParams = Fitting.sanitize(params(iSet));
        
    % TODO - smarter setting of seed?
    thisParams.seed = randi(1000000000);
    
    nTrials = length(choices{iSet});
    
    % Repeat signals to vectorize 'nInner' as if it were more trials
    signals{iSet} = repmat(signals{iSet}, nInner, 1);
    
    sim_results = Model.runVectorized(thisParams, signals{iSet});
    prob_choice = 1 ./ (1 + exp(-sim_results.lpo(:, end) / thisParams.temperature));
    
    % Reshape prob_choice from model to [nTrials x nInner]
    prob_choice = reshape(prob_choice, nTrials, nInner);
    
    % likelihood_samples is [nTrials x nInner] and contains the full choice likelihood for each run of
    % the model on each choice.
    chosePos = choices{iSet} == +1;
    choseNeg = ~chosePos;
    likelihood_samples{iSet}(chosePos, :) = thisParams.lapse/2 + (1-thisParams.lapse) * prob_choice(chosePos, :);
    likelihood_samples{iSet}(choseNeg, :) = thisParams.lapse/2 + (1-thisParams.lapse) * (1-prob_choice(choseNeg, :));
    
    log_likelihood_per_trial{iSet} = log(nanmean(likelihood_samples{iSet}, 2));
    log_likelihood = log_likelihood + sum(log_likelihood_per_trial{iSet});
end    

end