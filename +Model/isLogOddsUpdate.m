function [lpo, x_samples, weights] = isLogOddsUpdate(params, e, lpo)
%MODEL.ISLOGODDSUPDATE compute update to log odds of C using importance sampling model.

trials = size(e, 1);
oz = ones(trials, 1);

updates = params.updates;
noise = params.noise;
gamma = params.gamma;
sig_s = sqrt(params.var_s);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
samples = params.samples;

% Create two distributions representing p(x|C=+1) and p(x|C=-1) in the generative model
p_x_Cp = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]);
p_x_Cm = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]);

x_samples = zeros(trials, samples, updates);
weights = zeros(trials, samples, updates);

for n=1:updates
    % Convert from lpo (log odds) to the probability that C=+1
    pC = 1 ./ (1 + exp(-lpo));
    
    % Create likelihoods. Format is a matrix where each row specifies triples of [mu, sigma, pi] of
    % a mixture of Gaussians. Only one mode in the likelihood, but it's useful to use the MoG
    % format. See @mog.create
    likelihoods = [e(:) sig_s*oz oz];
    
    % Create the prior on x by marginalizing over the current posterior of C. The prior is also a
    % mixture of gaussians, but with 2 modes corresponding to C = +/-1
    p = p_match * pC + (1 - p_match) * (1 - pC);
    priors = [+oz, sig_x*oz, p, -oz, sig_x*oz, 1-p];
        
    % Q is the distribution from which samples of x are drawn; it is the current estimate of the
    % posterior over x using lpo as the prior over C
    Q = mog.prod(likelihoods, priors);
    
    % Draw samples from Q
    samp = mog.sample(Q, samples);
    x_samples(:, :, n) = samp;
    flat_samples = samp(:);
    
    % Get unnormalized importance weights for each sample using vectorized function mog.pdf
    weights(:, :, n) = 1 ./ mog.pdf(x_samples(:, :, n), priors);

    % Normalize importance-sampling weights
    if params.importance_norm
        weights(:, :, n) = weights(:, :, n) ./ sum(weights(:, :, n), 2);
    end

    % Compute p(x|C=+1) and p(x|C=-1), then take weighted sum for each trial.
    pCp = sum(reshape(mog.pdf(flat_samples, p_x_Cp), [trials samples]) .* weights(:, :, n), 2);
    pCm = sum(reshape(mog.pdf(flat_samples, p_x_Cm), [trials samples]) .* weights(:, :, n), 2);

    % Log likelihood odds is log(pCp/pCm)
    llo = log(pCp) - log(pCm);
    
    lpo = lpo * (1 - gamma / updates) + llo / updates;
    
    % Add zero-mean additive noise.
    lpo = lpo + randn(trials, 1) * noise;
end

end