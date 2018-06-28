function results = runLatentZNormalX(params, data)
% Variational Bayes inference in a model where C=+/-1, p(x|C) is a gaussian centered on +/-C with
% some probability, and p(s|x) is also gaussian. This version assumes a posterior that factorizes 
% x,z from c, i.e.
%    q(x,z,c)=q(x,z)*q(c)

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

pC = 1 ./ (1 + exp(-log_pC(:, 1)));
mu_x_trace(:, 1, 1) = 2 * pC - 1;
mu_x_trace(:, 1, 2) = 1 - 2 * pC;
log_odds_z_trace(:, 1) = log_prior_odds_z;

for t=1:frames
    % Implement VB to get t+1 probabilities
    log_odds_C = log_pC(:,t); % Use last posterior as new prior
    
    for i=1:updates
        pC = 1 ./ (1 + exp(-log_odds_C));
        mu_C = 2 * pC - 1; % C is bernoulli but with values {-1, +1}
        
        % The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
        mu_x_pos = (data(:,t) * var_x + mu_C * var_s) / var_xs; % z = +1
        mu_x_neg = (data(:,t) * var_x - mu_C * var_s) / var_xs; % z = -1

        % pz is the mass in each mode of the MOG
        log_odds_z = log_prior_odds_z + 2 * data(:, t) .* mu_C / var_xs;
        pz = 1 ./ (1 + exp(-log_odds_z));
        
        % Mean of x is weighted sum of means from each mode
        mu_x = mu_x_pos .* pz + mu_x_neg .* (1 - pz);
        
        % Infer updated p(c|s) and loop to next update for x,z now using log_odds_C for q(C)
        % TODO - for better consistency with sampling model, should we change this update to
        % log_odds_C on the RHS with gamma? Is 'updates' more like 'batch' or like 'samples'?
        log_odds_C = 2 * mu_x / var_x + log_pC(:,t);
    end
    
    % Result is 'new' log odds optionally with 'old' log odds subtracted out
    log_pC(:, t+1) = log_odds_C - gamma * log_pC(:,t);

    % Store 'trace' of x and z
    mu_x_trace(:, t+1, :) = [mu_x_pos(:) mu_x_neg(:)];
    log_odds_z_trace(:, t+1) = log_odds_z;
    
    % Add multiplicative noise to accumulated log probability
    if noise > 0
        log_pC(:, t+1) = log_pC(:, t+1) .* exp(noise * randn(trials,1) - noise^2/2);
    end
end

% compute decision variable
choices = log_pC(:, frames+1) > 0;

% Store results in the 'results' struct
results.walk = log_pC;
results.params = params;
results.choices = choices;
end
