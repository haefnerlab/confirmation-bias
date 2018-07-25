function results = runLatentZNormalX(params, data)
% Variational Bayes inference in a model where C=+/-1, p(x|C) is a gaussian centered on +/-C with
% some probability, and p(s|x) is also gaussian. This version assumes a posterior that factorizes 
% x,z from c, i.e.
%    q(x,z,c)=q(x,z|c)*q(c)
%...that is, q(x_t,z_t|c) is conditioned only on past values of c; there is no "lookahead" into the
% future as would be required for exact inference of q(x_t,z_t|data)

if ~exist('data', 'var')
    data = SamplingModel.genDataWithParams(params);
end

trials = params.trials;
frames = params.frames;
var_s = params.var_s;
var_x = params.var_x;
noise = params.noise;
gamma = params.gamma;
updates = params.updates;
pz0 = params.p_match;

log_prior_odds_z = log(pz0) - log(1 - pz0);
var_xs = var_x + var_s;

log_pC = zeros(trials, frames + 1);
log_pC(:, 1) = log(params.prior_C) - log(1 - params.prior_C);

mu_x_trace = zeros(trials, frames + 1, 2);
log_odds_z_trace = zeros(trials, frames + 1);

pC = sigmoid(log_pC(:, 1));
mu_x_trace(:, 1, 1) = bernoulli_plusminus(pC);
mu_x_trace(:, 1, 2) = bernoulli_plusminus(1-pC);
log_odds_z_trace(:, 1) = log_prior_odds_z;

for t=1:frames
    % Implement VB to get t+1 probabilities
    log_odds_C = log_pC(:,t); % Use last posterior as new prior
    
    % Analogously to the 'number of batches' option in the sampling model, here we assume that
    % q(x,z) and q(c) are potentially updated multiple times per stimulus frame. To implement this,
    % we take the term that would be the update for a single loop ('update_to_log_C' below) and
    % divide it by the number of updates. While this is no longer exact VB, it closely approximates
    % the exact solution and is more plausible for an online system that lacks the ability to
    % "buffer".
    for i=1:updates
        pC = sigmoid(log_odds_C);
        mu_C = bernoulli_plusminus(pC);
        
        % The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
        mu_x_pos = (data(:,t) * var_x + mu_C * var_s) / var_xs; % z = +1
        mu_x_neg = (data(:,t) * var_x - mu_C * var_s) / var_xs; % z = -1

        % pz is the mass in each mode of the MOG
        log_odds_z = log_prior_odds_z + 2 * data(:, t) .* mu_C / var_xs;
        pz = sigmoid(log_odds_z);
        
        % Mean of x is weighted sum of means from each mode
        mu_x = mu_x_pos .* pz + mu_x_neg .* (1 - pz);
        
        % Infer updated p(c|s) and loop to next update for x,z now using log_odds_C for q(C).
        % Note that log_pC is essentially a "baseline" value of log_odds_C; it stands in for
        % the prior term that would normally go there.
        update_to_log_C = 2 * mu_x / var_x;
        log_odds_C = log_odds_C * (1 - gamma / updates) + update_to_log_C / updates;
    end
    
    % Record posterior of C at frame t (recall index 1 is the prior, so index t+1 is for frame t)
    log_pC(:, t+1) = log_odds_C;

    % Store 'trace' of x and z
    mu_x_trace(:, t+1, :) = [mu_x_pos(:) mu_x_neg(:)];
    log_odds_z_trace(:, t+1) = log_odds_z;
    
    % Add multiplicative noise to accumulated log probability
    if noise > 0
        log_pC(:, t+1) = log_pC(:, t+1) .* exp(noise * randn(trials,1) - noise^2/2);
    end
end

% Compute decisions from sign of log posterior odds (optimal rule; no probability matching strategy)
choices = log_pC(:, frames+1) > 0;

% Store results in the 'results' struct
results.walk = log_pC;
results.params = params;
results.choices = choices;
end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function mu = bernoulli_plusminus(p)
% If x \in {-1, +1}, and p is the probability that x is +1, then E[x] is:
mu = 2 * p - 1;
end