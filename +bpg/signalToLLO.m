function [llo, emp_sensory_info] = signalToLLO(signal, kappa, p_match, add_noise_sig)
%SIGNALTOLLO take the output of bpg.getSignal and convert it to log-likelihood-odds using a
%precomputed histogram of possible signal values per noise level. The second argument (kappa)
%controls behavior: if kappa is empty, it is treated as 'unknown' and the signal histogram is
%treated as a mixture of distributions (modeling 'uknown noise level'). If kappa is a scalar, only
%the histogram for that signal level is used.
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
end

if isempty(kappa)
    % The underlying signal level is uknown. Assume a flat prior over kappas in the above table.
    likelihood_positive_per_mode = normpdf(signal(:)', +kappa_mu_sig(:,2), kappa_mu_sig(:,3));
    likelihood_negative_per_mode = normpdf(signal(:)', -kappa_mu_sig(:,2), kappa_mu_sig(:,3));
    
    likelihood_positive = p_match * mean(likelihood_positive_per_mode,1) + (1 - p_match) * mean(likelihood_negative_per_mode,1);
    likelihood_negative = p_match * mean(likelihood_negative_per_mode,1) + (1 - p_match) * mean(likelihood_positive_per_mode,1);
    
    sensory_llo = log(mean(likelihood_positive_per_mode,1)) - log(mean(likelihood_negative_per_mode,1));
else
    % For each given signal value, look up parameters of the gaussian fit for its corresponding signal
    % level (kappa)
    signal_mu = interp1(kappa_mu_sig(:,1), kappa_mu_sig(:,2), kappa);
    signal_sig = interp1(kappa_mu_sig(:,1), kappa_mu_sig(:,3), kappa);

    % For mixtures of Gaussians, log-likelihood odds are not a nice function of the log likelihood of
    % the corresponding Gaussian components. There is no real shortcut besides computing the likelihood
    % directly then taking the log.
    likelihood_positive = p_match * normpdf(signal(:)', +signal_mu, signal_sig) + (1 - p_match) * normpdf(signal(:)', -signal_mu, signal_sig);
    likelihood_negative = p_match * normpdf(signal(:)', -signal_mu, signal_sig) + (1 - p_match) * normpdf(signal(:)', +signal_mu, signal_sig);
    
    sensory_llo = log(normpdf(signal(:)', +signal_mu, signal_sig)) - log(normpdf(signal(:)', -signal_mu, signal_sig));
end

llo = log(likelihood_positive) - log(likelihood_negative);
emp_sensory_info = 1 ./ (1 + exp(-abs(sensory_llo)));

llo = reshape(llo, size(signal));
emp_sensory_info = reshape(emp_sensory_info, size(signal));

end