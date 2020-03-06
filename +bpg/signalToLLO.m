function [full_llo, sensory_llo, emp_sensory_info] = signalToLLO(signal, kappa, p_match, add_noise_sig)
%SIGNALTOLLO take the output of bpg.getSignal and convert it to log-likelihood-odds using a
%precomputed histogram of possible signal values per noise level. The second argument (kappa)
%controls behavior: if kappa is empty, it is treated as 'unknown' and the signal histogram is
%treated as a mixture of distributions (modeling 'uknown noise level'). If kappa is a scalar, only
%the histogram for that signal level is used. Finally, both 'kappa' and 'p_match' may be [M x 2]
%vectors indicating discrete distributions of values, i.e. the subject's 'prior' over each unknown
%value. For instance kappa = [0 .6; 0.04 .4] means kappa has a 60% (resp. 40%) chance of having a
%value of 0 (resp. 0.04);
%
%Also estimates "empirical sensory info" as the likelihood ratio if p_match had been equal to 1
%(i.e. with all uncertainty reduced to the sensory kind). The relevant quantity for the model is,
%instead, "expected sensory info" over many trials.

% We precomputed the relationship between stimulus noise (kappa) and empirical distributions
% of signal values, as returned by 'bpg.getSignal'. Using the table directly here:

kappa_mu_sig = [
    0.0000    0.0000    0.1045
    0.0400    0.0864    0.1047
    0.0800    0.3521    0.2138
    0.1200    0.7818    0.3263
    0.1600    1.3820    0.4471
    0.2000    2.1585    0.5881
    0.2400    3.0953    0.7398
    0.2800    4.2235    0.9350
    0.3200    5.4120    1.1071
    0.3600    6.8388    1.3741
    0.4000    8.3841    1.5789
    0.4400   10.2649    1.9241
    0.4800   11.4003    2.0046
    0.5200   13.6046    2.3481
    0.5600   15.7614    2.7215
    0.6000   17.3269    2.9511
    0.6400   18.4907    3.0103
    0.6800   21.2790    3.4385
    0.7200   22.6857    3.6132
    0.7600   24.3695    4.0899
    0.8000   26.2707    4.0179];

if exist('add_noise_sig', 'var')
    % Additional sensory noise is equated with added variance to the signal distribution
    kappa_mu_sig(:,3) = sqrt(kappa_mu_sig(:,3).^2 + add_noise_sig^2);
    signal = signal + add_noise_sig * randn(size(signal));
end

%% Convert input arguments in whatever form into a prior over 'kappa' and 'p_match'
if isempty(kappa)
    kappa_prior = [kappa_mu_sig(:,1), ones(size(kappa_mu_sig(:,1)))];
elseif isscalar(kappa)
    kappa_prior = [kappa 1];
else
    kappa_prior = kappa;
end

if isempty(p_match)
    values = 0.5:0.1:1;
    p_match_prior = [values(:) ones(size(values(:)))];
elseif isscalar(p_match)
    p_match_prior = [p_match 1];
else
    p_match_prior = p_match;
end

% Ensure normalization
kappa_prior(:,2) = kappa_prior(:,2)/sum(kappa_prior(:,2));
p_match_prior(:,2) = p_match_prior(:,2)/sum(p_match_prior(:,2));

% Interpolate the above table to get mu and sigma parameters of signal distributions at points
% specified in the prior over kappa
prior_mu_per_kappa = interp1(kappa_mu_sig(:,1), kappa_mu_sig(:,2), kappa_prior(:,1));
prior_sig_per_kappa = interp1(kappa_mu_sig(:,1), kappa_mu_sig(:,3), kappa_prior(:,1));

K = size(kappa_prior, 1);
P = size(p_match_prior, 1);
S = numel(signal);

%% Compute a likelihood term per-kappa and per-p_match (independent priors)

% Each of these terms is [K x S] where K is number of 'kappa' values and 'S' is numel(signal)
likelihood_positive_per_kappa = normpdf(signal(:)', +prior_mu_per_kappa, prior_sig_per_kappa);
likelihood_negative_per_kappa = normpdf(signal(:)', -prior_mu_per_kappa, prior_sig_per_kappa);

% Compute 'sensory_llo' as log likelihood odds with p_match = 1
sensory_llo = log(likelihood_positive_per_kappa' * kappa_prior(:,2)) - ...
    log(likelihood_negative_per_kappa' * kappa_prior(:,2));
emp_sensory_info = 1 ./ (1 + exp(-abs(sensory_llo)));

% Expand to shape [P x K x S] by 'mixing' using p_match terms
likelihood_positive = p_match_prior(:,1) .* reshape(likelihood_positive_per_kappa, [1 K S]) + ...
    (1 - p_match_prior(:,1)) .* reshape(likelihood_negative_per_kappa, [1 K S]);
likelihood_negative = p_match_prior(:,1) .* reshape(likelihood_negative_per_kappa, [1 K S]) + ...
    (1 - p_match_prior(:,1)) .* reshape(likelihood_positive_per_kappa, [1 K S]);

% Marginalize over 2 prior dimensions to get [S x 1] marginal likelihoods per category
pr_k = reshape(kappa_prior(:,2), [1 K 1]);
pr_p = reshape(p_match_prior(:,2), [P 1 1]);
marginal_likelihood_positive = squeeze(sum(sum(likelihood_positive .* pr_k .* pr_p, 1), 2));
marginal_likelihood_negative = squeeze(sum(sum(likelihood_negative .* pr_k .* pr_p, 1), 2));

full_llo = log(marginal_likelihood_positive) - log(marginal_likelihood_negative);

% Reshape output to match shape of 'signal'
full_llo = reshape(full_llo, size(signal));
sensory_llo = reshape(sensory_llo, size(signal));
emp_sensory_info = reshape(emp_sensory_info, size(signal));

end