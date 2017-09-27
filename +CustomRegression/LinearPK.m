function [sob, postVal, errors] = LinearPK(data, responses, standardize)
%LinearPK Regress PK as a linear function.
%
% sob = LinearPK(data, responses) returns sob = [slope offset bias].

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

    function NLL = neg_bernoulli_log_likelihood(sob)
        weights = sob(2) + (0:frames-1) * sob(1);
        logits = data * weights(:) + sob(3);
        neg_log_bernoulli = -responses .* logits(:) + log(1 + exp(logits(:)));
        NLL = sum(neg_log_bernoulli);
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
