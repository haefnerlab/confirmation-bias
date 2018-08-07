function lpo = vbLogOddsUpdate(params, e, lpo)
%MODEL.VBLOGODDSUPDATE compute log likelihood odds estimate using variational Bayes model,
%vectorized over trials. Uses q(x,z|C)q(C) factorization.

trials = size(e, 1);

updates = params.updates;
noise = params.noise;
gamma = params.gamma;
var_s = params.var_s;
var_x = params.var_x;
var_xs = var_s + var_x;

pz0 = params.p_match;
log_prior_odds_z = log(pz0) - log(1 - pz0);

% Initialize q(C) to the running posterior
q_odds_C = lpo;

for n=1:updates
    % Convert from lpo (log odds) to the expected value of C
    mu_C = bernoulli_plusminus(sigmoid(q_odds_C));
    
    % The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
    mu_x_pos = (e * var_x + mu_C * var_s) / var_xs; % z = +1
    mu_x_neg = (e * var_x - mu_C * var_s) / var_xs; % z = -1
    
    % pi_z is the mass in the z=+1 mode of the MoG
    log_odds_z = log_prior_odds_z + 2 * e .* mu_C / var_xs;
    pi_z = sigmoid(log_odds_z);
    
    % Compute updated log odds of C
    q_odds_C = lpo + 2 * (pi_z .* mu_x_pos - (1 - pi_z) .* mu_x_neg) / var_x;
    
    % Add noise (multiply by log-normal random variable with expected value 1)
    eta = exp(randn(trials, 1) * noise - noise^2/2);
    q_odds_C = eta .* q_odds_C;
end

lpo = q_odds_C - gamma * lpo;

end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function mu = bernoulli_plusminus(p)
% If x \in {-1, +1}, and p is the probability that x is +1, then E[x] is:
mu = 2 * p - 1;
end