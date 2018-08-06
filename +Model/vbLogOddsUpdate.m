function [llo] = vbLogLikelihood(params, e, lpo)
%MODEL.VBLOGLIKELIHOOD compute log likelihood odds estimate using variational Bayes model,
%vectorized over trials. Uses q(x,z|C)q(C) factorization

var_s = params.var_s;
var_x = params.var_x;
var_xs = var_s + var_x;

pz0 = params.p_match;
log_prior_odds_z = log(pz0) - log(1 - pz0);

% Convert from lpo (log odds) to the probability that C=+1
pC = 1 ./ (1 + exp(-lpo));
mu_C = bernoulli_plusminus(pC);
        
% The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
mu_x_pos = (e * var_x + mu_C * var_s) / var_xs; % z = +1
mu_x_neg = (e * var_x - mu_C * var_s) / var_xs; % z = -1

% pz is the mass in each mode of the MOG
log_odds_z = log_prior_odds_z + 2 * e .* mu_C / var_xs;
pz = sigmoid(log_odds_z);

% Mean of x is weighted sum of means from each mode
mu_x = mu_x_pos .* pz + mu_x_neg .* (1 - pz);

% Infer updated p(c|s) and loop to next update for x,z now using log_odds_C for q(C). Note that the
% full VB update is derived as log(p(C)) <-- log(p(C)) + 2*mu_x/var_x . We do not include the
% log(p(C)) term on the RHS since that is handled by the outer loop.
llo = 2 * mu_x / var_x;

end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function mu = bernoulli_plusminus(p)
% If x \in {-1, +1}, and p is the probability that x is +1, then E[x] is:
mu = 2 * p - 1;
end