function [log_likelihood, log_likelihood_per_trial, likelihood_samples] = choiceModelLogLikelihood(params, signals, choices, nInner)

params.var_s = params.var_s_per_sample * params.samples;

% TODO - smarter setting of seed?
params.seed = randi(1000000000);

nTrials = length(choices);

% Repeat signals to vectorize 'nInner' as if it were more trials
signals = repmat(signals, nInner, 1);

sim_results = Model.runVectorized(params, signals);
prob_choice = 1 ./ (1 + exp(-params.alpha * sim_results.lpo(:, end)));

% Reshape prob_choice from model to [nTrials x nInner]
prob_choice = reshape(prob_choice, nTrials, nInner);

% likelihood_samples is [nTrials x nInner] and contains the full choice likelihood for each run of
% the model on each choice.
chosePos = choices == +1;
choseNeg = ~chosePos;
likelihood_samples(chosePos, :) = params.lapse/2 + (1-params.lapse) * prob_choice(chosPos, :);
likelihood_samples(choseNeg, :) = params.lapse/2 + (1-params.lapse) * (1-prob_choice(choseNeg, :));

log_likelihood_per_trial = log(mean(likelihood_samples, 2));
log_likelihood = sum(log_likelihood_per_trial);

end