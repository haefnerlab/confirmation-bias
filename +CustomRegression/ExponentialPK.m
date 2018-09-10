function [abb, postVal, errors] = ExponentialPK(data, responses, standardize)
%EXPONENTIALPK Regress PK as an exponential function of time.
%
% abb = LinearPK(data, responses) returns abb = [alpha beta bias] such that w=alpha*exp(beta*f)
% and choices~sigmoid(signal'*w+bias)

if nargin < 3, standardize = 0; end

% Standardize each regressor.
switch standardize
    case 0
        % do nothing
    case 1
        % assume 0 mean (nothing to subtact) and iid (std taken over all data)
        data = data / std(data(:));
    case 2
        data = zscore(data);
    otherwise
        error('Expected argument ''standardize'' to be one of [0, 1, 2]');
end

% convert boolean to float type
assert(islogical(responses));
responses = 1.0 * responses(:);

[~, frames] = size(data);

    function NLL = neg_bernoulli_log_likelihood(abb)
        weights = abb(1) * exp(abb(2) * (0:frames-1));
        logits = data * weights(:) + abb(3);
        neg_log_bernoulli = -responses .* logits(:) + log(1 + exp(logits(:)));
        NLL = sum(neg_log_bernoulli);
        NLL = NLL + 1/2*abb(2)^2/100;
    end

compute_error = nargout > 2;

glm_weights = glmfit(data, responses, 'binomial');
ab = CustomRegression.expFit(glm_weights(2:end));
init_guess = [ab glm_weights(1)];

% Fit weights using 'fminunc', only computing the hessian (which is slow) if errors are requested.
if compute_error
    [abb, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_bernoulli_log_likelihood, init_guess);
    % attempt to invert the hessian for standard error estimate - this sometimes fails silently,
    % returning NaN.
    errors = sqrt(diag(abs(inv(-hessian))));
else
    [abb, negPostVal] = fminunc(@neg_bernoulli_log_likelihood, init_guess);
end
postVal = -negPostVal;
end
