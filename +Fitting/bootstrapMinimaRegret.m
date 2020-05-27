function [n99, pr_n] = bootstrapMinimaRegret(minima, stdevs, tolerance)
%FITTING.BOOTSRAPMINIMAREGRET bootstrap the given set of minima to estimate the probability that
%they contain a global minimum. Use this to estimate 'n99', the minimum # searches required to be
%99% sure that one of them contains the global minimum. Intuitively, if most of the given minima
%agree, then it's likely they contain the global minimum. If, say, one of the minima is far better
%than the others, then chances are the space has not yet been explored well. This case manifests as
%a high probability of 'regret' after resampling the minima with replacement.
%
% [n99, pr_n] = FITTING.BOOTSRAPMINIMAREGRET(minima) given a set of n local minima from some search
% procedure, computes 'pr_n', a [length(minima) x 1] vector containing 1 minus the  "regret"
% probability for each value of n from 1:length(minima). Regret is defined as being farther than
% 'tolerance' away from the minimum of all the given minima, where default tolerance is 1. pr_n is
% then extrapolated with an exponential fit to esimate where it crosses 99%. That is, n99 is the
% estimated number of searches it would take to be 99% sure that there is no better minimum
% somewhere else out there.
%
% [n99, pr_n] = FITTING.BOOTSRAPMINIMAREGRET(minima, stdevs) also specify uncertainty around each of
% the given minima. The bootstrapping procedure then adds 'stdev' amount of noise on each run before
% computing regret. In this case, it may be sensible to lower the tolerance to 0. By default (or if
% set to []), stdevs is 0.
%
% [n99, pr_n] = FITTING.BOOTSRAPMINIMAREGRET(minima, stdevs, tolerance) specifies how close is
% close-enough to be unregretful. For instance, if min(minima) is 10 and tolerance is 1, then some
% other local minimum at 10.9 is considered close enough and does not count towards regret. Default
% tolerance is 1.

if nargin < 2 || isempty(stdevs), stdevs = zeros(size(minima)); end
if nargin < 3, tolerance = 1; end
if isscalar(stdevs), stdevs = repmat(stdevs, size(minima)); end

global_min = min(minima);

if all(stdevs == 0) && all(minima < global_min+tolerance)
    n99 = 1;
    pr_n = ones(size(minima));
    return;
end

for n=length(minima):-1:1
    for b=1000:-1:1
        boot_idx = randi(length(minima), 1, n);
        boot_minima = minima(boot_idx) + randn(1,n).*stdevs(boot_idx);
        % Get estimated 'regret', equal to the distance we are from the "true" global minimum under
        % the hypothetical where only n runs were collected. Always positive, unless 'stdevs' is
        % given, in which case this may be negative (accounting for probability that the _value_ of
        % each estimated minimum is lower than reported).
        boot_regret(n,b) = min(boot_minima) - global_min;
    end
end

% What is the probability (fraction of bootstraps) of having regret < tolerance? In other words,
% what's the probability with 'n' simulations that the best min we find is within 'tolerance' of the
% "global" min? Higher pr_n means n is more likely to be sufficient to find "the" global min.
pr_n = mean(boot_regret < tolerance, 2);

% The expected shape of pr_n is an exponential CDF. Get a quick estimate of the mean parameter and
% use it to extrapolate (or interpolate) the lowest 'n' at which pr(regret < tolerance) > 99% . Sum
% of squares is not exactly the right noise model, but whatever.
exp_mu = fminsearch(@(mu) sum((1-exp(-(1:length(minima))/mu)-pr_n').^2), length(minima)/2);

% Solve for where exp cdf >= 0.99
n99 = ceil(expinv(0.99, exp_mu));
end