function [weights, postVal, errors] = PsychophysicalKernel(data, responses, hpr_ridge, hpr_ar1, split_smoothness)
%PSYCHOPHYSICALKERNEL Regress PK
%
% [ weights, errors ] = PSYCHOPHYSICALKERNEL(data, responses) where data is
% [trials x regressors] and responses is [trials x 1] of booleans.

if nargin < 5, split_smoothness = false; end
if nargin < 4, hpr_ar1 = 0; end
if nargin < 3, hpr_ridge = 0; end

% convert boolean to float type
responses = 1.0 * responses;

[~, p] = size(data);

break_points = p-1; % don't smooth onto bias term
if split_smoothness, break_points = [break_points floor(p/2)]; end
D = derivative_matrix(p, break_points);

negLogPost = @(w) -log_prior(w, D, hpr_ridge, hpr_ar1) - bernoulli_log_likelihood(data, responses, w);

if nargout > 2
    [weights, negPostVal, ~, ~, ~, hessian] = fminunc(negLogPost, zeros(p, 1));
    % attempt to invert the hessian for error bars
    errors = abs(diag(inv(hessian)));
else
    [weights, negPostVal] = fminunc(negLogPost, zeros(p, 1));
end

postVal = -negPostVal;

end

function D = derivative_matrix(n, break_at)
D = - eye(n) + diag(ones(n-1, 1), 1);
D(end) = 0;
for i=break_at
    D(i,i) = 0;
    D(i, i+1) = 0;
end
end

function LP = log_prior(weights, D, hpr_ridge, hpr_ar1)
ridge = -0.5 * dot(weights, weights);
ar1 = 0;
if hpr_ar1 > 0
    ar1 = -0.5 * dot(D*weights, D*weights);
end
LP = hpr_ridge * ridge + hpr_ar1 * ar1;
end

function LL = bernoulli_log_likelihood(data, responses, weights)
p = sigmoid(data * weights);
LL = dot(responses, log(p)) + dot(1-responses, log(1-p));
end

function y = sigmoid(x)
y = (1 + exp(-x)).^-1;
end