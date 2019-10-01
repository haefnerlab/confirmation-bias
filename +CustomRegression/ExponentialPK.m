function [abb, postVal, errors] = ExponentialPK(data, responses, standardize)
%EXPONENTIALPK Regress PK as an exponential function of time.
%
% abb = EXPONENTIALPK(data, responses) returns abb = [alpha beta bias] such that w=alpha*exp(beta*f)
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

    function [NLP, gradNLP] = neg_log_posterior(abb)
        weights = abb(1) * exp(abb(2) * (0:frames-1));
        logits = data * weights(:) + abb(3);
        neg_log_bernoulli = -responses .* logits(:) + log(1 + exp(logits(:)));
        NLP = sum(neg_log_bernoulli);
        % Add mild priors to alpha and beta terms for stability
        neg_log_prior_beta = 1/2*abb(2)^2/100;
        neg_log_prior_alpha = abs(abb(1))/100;
        NLP = NLP + neg_log_prior_alpha + neg_log_prior_beta;
        
        if nargout >= 2
            dNLP_dlogits = -responses + exp(logits(:))./(1+exp(logits(:)));
            dlogits_dw = data';
            dw_dalpha = exp(abb(2)*(0:frames-1));
            dw_dbeta = abb(1)*(0:frames-1).*exp(abb(2)*(0:frames-1));
            gradNLP = [dw_dalpha * dlogits_dw * dNLP_dlogits + 1/10, ...
                dw_dbeta * dlogits_dw * dNLP_dlogits + abb(2)/10, ...
                sum(dNLP_dlogits)];
        end
    end

compute_error = nargout > 2;

init_guess = [1 0 0];

% Fit weights using 'fminunc', only computing the hessian (which is slow) if errors are requested.
options = optimoptions('fminunc', ...
    'Display', 'none', ...
    'Algorithm', 'trust-region', ...
    'SpecifyObjectiveGradient', true);

if compute_error
    [abb, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_log_posterior, init_guess, options);
    % attempt to invert the hessian for standard error estimate - this sometimes fails silently,
    % returning NaN.
    errors = sqrt(diag(abs(inv(-hessian))));
else
    [abb, negPostVal] = fminunc(@neg_log_posterior, init_guess, options);
end
postVal = -negPostVal;
end
