function results = runLatentZNormalX(params)
% Variational Bayes inference in a model where C=+/-1, p(x|C) is a gaussian centered on +/-C with
% some probability, and p(s|x) is also gaussian. This version assumes a posterior that factorizes 
% x,z from c, i.e.
%    q(x,z,c)=q(x,z)*q(c)

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
    prior_C = pC(:,t); % Use last posterior as new prior
    log_prior_odds_C = log(prior_C) - log(1 - prior_C);
    
    for i=1:updates
        % The form of q(x,z) is a mixture of two gaussians corresponding to z = +/-1
        mu_pos = (data(:,t) * var_x + (2 * prior_C - 1) * var_s) / (var_s + var_x); % z = +1
        mu_neg = (data(:,t) * var_x - (2 * prior_C - 1) * var_s) / (var_s + var_x); % z = -1

        % pz is the mass in each mode of the MOG
        log_odds_z = log(pz0) - log(1-pz0) + 2 * data(:, t) .* (2 * prior_C - 1);
        pz = 1 ./ (1 + exp(-log_odds_z));
        
        % Mean of x is weighted sum of means from each mode
        mu_x = mu_pos .* pz + mu_neg .* (1 - pz);
        
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
results.log_odds = pC;
results.params = params;
results.choices = choices;
end
