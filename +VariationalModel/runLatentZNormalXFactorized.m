function results = runLatentZNormalXFactorized(params)
% Variational Bayes inference in a model where C=+/-1, p(x|C) is a gaussian centered on +/-C with
% some probability, and p(s|x) is also gaussian. This version assumes a posterior that factorizes
% x,z and c, i.e.
%    q(x,z,c)=q(x)*q(z)*q(c)

trials = params.trials;
frames = params.frames;
var_s = params.var_s;
var_x = params.var_x;
noise = params.noise;
gamma = params.gamma;
updates = params.updates;
pz0 = params.p_match;

pC = zeros(trials, frames + 1);
pC(:, 1) = params.prior_C;

% Generate the stimulus, all with C=+1 'correct'
data = SamplingModel.genDataWithParams(params);

for t=1:frames
    % Implement VB to get t+1 probabilities
    prior_C = pC(:, t); % Use last posterior as new prior
    log_prior_odds_C = log(prior_C) - log(1 - prior_C);
    pz = repmat(pz0, trials, 1);
    
    for i=1:updates
        % Infer p(x|z,c,s)
        expected_cz = (2 * prior_C - 1) .* (2 * pz - 1);
        mu_x = (expected_cz * var_s + data(:, t) * var_x) / (var_s + var_x);
        % Note: updated var_x is constant: (1/var_x + 1/var_s)^-1

        % Infer p(z|x,c)
        log_odds_z = 2 * (2 * prior_C - 1) .* mu_x / var_x + log(pz0) - log(1-pz0);
        pz = 1 ./ (1 + exp(-log_odds_z));
        
        % Infer updated p(c|s)
        log_odds_C = 2 * (2 * pz - 1) .* mu_x / var_x + log_prior_odds_C * (1 - gamma);
        
        % Add multiplicative noise to accumulated log probability
        if noise > 0
            log_odds_C = log_odds_C .* exp(noise * randn(trials,1));
        end
        pC(:, t+1) = 1 ./ (1 + exp(-log_odds_C)); 
    end
end

% compute decision variable
choices = pC(:, frames+1) > .5;

% Store results in the 'results' struct
results.walk = pC;
results.params = params;
results.choices = choices;
end