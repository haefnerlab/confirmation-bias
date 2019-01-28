function lpo = vbLogOddsUpdateCZX(params, e, lpo)
%MODEL.VBLOGODDSUPDATECZX compute log likelihood odds estimate using variational Bayes model,
%vectorized over trials. Uses full q(x)q(C)q(z) factorization.

trials = size(e, 1);

updates = params.updates;
noise = params.noise;
gamma = params.gamma;
var_s = params.var_s;
step_size = params.step_size;
var_x = params.var_x;
var_xs = var_s + var_x;

pz0 = params.p_match;
log_prior_odds_z = log(pz0) - log(1 - pz0);

% Initialize mu_z to the prior
mu_z = bernoulli_plusminus(sigmoid(log_prior_odds_z));

for n=1:updates
    % Convert from lpo (log odds) to the expected value of C
    mu_C = bernoulli_plusminus(sigmoid(lpo));
    
    % Update x
    mu_x = (var_s * mu_C .* mu_z + var_x * e) / var_xs;
    
    % Update z
    log_odds_z = log_prior_odds_z + 2 * mu_x .* mu_C / var_xs;
    mu_z = bernoulli_plusminus(sigmoid(log_odds_z));
    
    % Update C using (gamma-discounted) prior and (step_size-discounted) update rule
    lpo = lpo * (1 - gamma / updates) + step_size * 2 * mu_x .* mu_z / var_x / updates;
    
    % Add zero-mean additive noise.
    lpo = lpo + randn(trials, 1) * noise;
end

end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function mu = bernoulli_plusminus(p)
% If x \in {-1, +1}, and p is the probability that x is +1, then E[x] is:
mu = 2 * p - 1;
end