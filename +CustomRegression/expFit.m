function [expfit, errors] = expFit(weights, errors)
if nargin < 2
    variances = ones(size(weights));
else
    variances = errors.^2;
end

% parameters are [scale beta] where the fit is weights = scale*exp(frames*beta). Negative beta is
% primacy, positive is recency.
p0 = [mean(weights) 0];

options = optimoptions(@fmincon, 'SpecifyObjectiveGradient', true);
[expfit, ~, ~, ~, ~, ~, hessian] = fmincon(@(p) errfn(p, weights(:), variances(:)), ...
    p0, [], [], [], [], [0 -10], [inf 10], [], options);
errors = sqrt(diag(abs(inv(-hessian))));

end

function [err, grad] = errfn(p, w, v)
scale = p(1);
beta = p(2);
betaPriorVar = 100;
frames = (0:length(w)-1)';
expvals = scale*exp(frames*beta);
err = sum((w - expvals).^2 ./ v) + beta^2 / betaPriorVar;

if nargout >= 2
    gradScale = 2 * sum((scale*exp(2*frames*beta) - w.*exp(frames*beta)) ./ v);
    gradBeta = 2 * sum(frames .* (expvals .* (expvals - w)) ./ v) + 2 * beta / betaPriorVar;
    grad = [gradScale gradBeta];
    
    % if nargout >= 3
    %     hessScale = 2 * sum(exp(2*frames*beta) ./ v);
    %     hessScaleBeta = 2 * sum((2*scale*exp(2*frames*beta) - w.*exp(frames*beta)) .* frames ./ v);
    %     hessBeta = 2 * sum((frames.^2 .* expvals .* (w + 2)) ./ v);
    %     hess = [hessScale hessScaleBeta; hessScaleBeta hessBeta];
    % end
end
end