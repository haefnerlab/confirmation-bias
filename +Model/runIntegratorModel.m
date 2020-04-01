function results = runIntegratorModel(params, signals)

[trials, frames] = size(signals);
results.params = params;

% By adjusting for temperature before integrating, we eliminate the dependency between the
% fit temperature and the fit prior + bound; signals are directly interpreted as log likelihood
% odds.
signals = signals / params.temperature;

% Initialize integrator to the log prior odds.
results.lpo = zeros(trials, frames+1);
results.lpo(:,1) = log(params.prior_C) - log(1-params.prior_C);

for f=1:frames
    results.lpo(:,f+1) = results.lpo(:,f)*(1-params.gamma) + signals(:, f);
    
    % Implement sticky bound by setting results.lpo to inf in the loop (it can never undo the
    % inf) and setting all inf values to the bound outside the loop.
    results.lpo(results.lpo >= +params.bound) = +inf;
    results.lpo(results.lpo <= -params.bound) = -inf;
end
results.lpo(isinf(results.lpo)) = sign(results.lpo(isinf(results.lpo)))*params.bound;

lapse_range = 1-(params.lapse_1+params.lapse_2);
results.prob_choice = params.lapse_1 + lapse_range * 1./(1+exp(-results.lpo(:,end)));
results.choice = rand(trials, 1) < results.prob_choice;
end