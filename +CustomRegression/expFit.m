function [expfit, errors] = expFit(weights, errors)
if nargin < 2
    variances = ones(size(weights));
else
    variances = errors.^2;
end

% parameters are [scale beta] where the fit is weights = scale*exp(frames*beta). Negative beta is
% primacy, positive is recency.
p0 = [mean(weights) 0];

options = optimoptions(@fmincon, 'SpecifyObjectiveGradient', true, 'Display', 'none');
[expfit, ~, ~, ~, ~, ~, hessian] = fmincon(@(p) errfn(p, weights(:), variances(:)), ...
    p0, [], [], [], [], [0 -10], [inf 10], [], options);
errors = sqrt(diag(inv(hessian)));

end

function [err, grad] = errfn(p, w, v)
scale = p(1);
beta = p(2);
betaPriorVar = 100;
frames = (0:length(w)-1)';
w_hat = scale*exp(frames*beta);
diffs = w_hat - w;
err = 1/2 * sum(diffs.^2 ./ v) + 1/2 * beta^2 / betaPriorVar;

if nargout >= 2
    gradScale = sum((scale*exp(2*frames*beta) - w.*exp(frames*beta)) ./ v);
    gradBeta = sum(frames .* (w_hat .* diffs) ./ v) + beta / betaPriorVar;
    grad = [gradScale gradBeta];
end
end