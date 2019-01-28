function lpo = vbLogOddsUpdate(params, e, lpo)
%MODEL.VBLOGODDSUPDATE compute log likelihood odds estimate using variational Bayes model,
%vectorized over trials. Uses q(x,z|C)q(C) factorization.

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

for n=1:updates
    % Convert from lpo (log odds) to the expected value of C
    mu_C = bernoulli_plusminus(sigmoid(lpo));
    
    % The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
    mu_x_pos = (e * var_x + mu_C * var_s) / var_xs; % z = +1
    mu_x_neg = (e * var_x - mu_C * var_s) / var_xs; % z = -1
    
    % pi_z is the mass in the z=+1 mode of the MoG
    log_odds_z = log_prior_odds_z + 2 * e .* mu_C / var_xs;
    pi_z = sigmoid(log_odds_z);
    
    % Compute updated log odds of C using (gamma-discounted) prior and (step_size-discounted) update
    % rule based on q(x,z)
    lpo = lpo * (1 - gamma / updates) + step_size * 2 * (pi_z .* mu_x_pos - (1 - pi_z) .* mu_x_neg) / var_x / updates;
    
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