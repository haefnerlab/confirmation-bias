function [sob, postVal, errors] = LinearPK(data, responses, data_mean)
%LinearPK Regress PK as a linear function.
%
% sob = LinearPK(data, responses) returns sob = [slope offset bias].

% Standardize each regressor.
if nargin < 3
    data = zscore(data);
else
    % We can do our own (more precise) z-scoring when we know the true
    % mean.
    data_centered = data - data_mean;
    % Note we can /N rather than /(N-1) for an unbiased variance estimate
    % here.
    data_variance = sum(data_centered.^2, 1) / size(data, 1);
    data = (data - data_mean) ./ sqrt(data_variance);
end

% convert boolean to float type
responses = 1.0 * responses;

[~, frames] = size(data);

    function NLL = neg_bernoulli_log_likelihood(sob)
        weights = sob(2) + (0:frames-1) * sob(1);
        p = sigmoid(data * weights' + sob(3))';
        NLL = -dot(responses, log(p)) - dot(1-responses, log(1-p));
    end

compute_error = nargout > 2;

% Fit weights using 'fminunc', only computing the hessian (which is
% slow) if errors are requested.
if compute_error
    [sob, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_bernoulli_log_likelihood, zeros(1, 3));
    % attempt to invert the hessian for standard error estimate - this
    % sometimes fails silently, returning NaN.
    errors = sqrt(diag(abs(inv(-hessian))));
else
    [sob, negPostVal] = fminunc(@neg_bernoulli_log_likelihood, zeros(1, 3));
end
postVal = -negPostVal;
end

function y = sigmoid(x)
y = (1 + exp(-x)).^-1;
end