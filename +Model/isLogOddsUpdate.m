function [llo, x_samples, weights] = isLogLikelihood(params, e, lpo)
%MODEL.ISLOGLIKELIHOOD compute log likelihood odds estimate using importance sampling model,
%vectorized over trials.

trials = size(e, 1);

sig_s = sqrt(params.var_s);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
samples = params.samples;

% Convert from lpo (log odds) to the probability that C=+1
pC = 1 ./ (1 + exp(-lpo));

for t=trials:-1:1
    % Likelihood is a unimodal Gaussian centered on e(t), but it is convenient to create a
    % mixture-of-Gaussians object to represent it.
    likelihood = mog.create(e(t), sig_s, 1);

    % Create the prior on x by marginalizing over the current posterior of C
    p = p_match * pC(t) + (1 - p_match) * (1 - pC(t));
    prior = mog.create([+1 -1], [sig_x sig_x], [p 1-p]);
    
    % Q is the distribution from which samples of x are drawn; it is the current estimate of the
    % posterior over x using lpo as the prior over C
    Q = mog.prod(likelihood, prior);
    
    % Draw samples from Q
    x_samples(t, :) = mog.sample(Q, samples);
    
    % Get unnormalized importance weights for each sample
    weights(t, :) = 1 ./ mog.pdf(x_samples(t, :), prior);
end

% Normalize importance-sampling weights
weights = weights ./ sum(weights, 2);

% Create two distributions representing p(x|C=+1) and p(x|C=-1) in the generative model
p_x_Cp = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]);
p_x_Cm = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]);

% Compute p(x|C=+1) and p(x|C=-1), then take weighted sum for each trial.
pCp = sum(mog.pdf(x_samples, p_x_Cp) .* weights, 2);
pCm = sum(mog.pdf(x_samples, p_x_Cm) .* weights, 2);

% Log likelihood odds is log(pCp/pCm)
llo = log(pCp) - log(pCm);

end