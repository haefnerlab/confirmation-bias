function expfit = expFit(weights, errors)
if nargin < 2, errors = ones(size(weights)); end

p0 = [mean(weights) 0 10];

expfit = fmincon(@(p) errfn(p, weights, errors), p0, [], [], [], [], [0 0 0], [inf inf inf]);

end

function [err, grad] = errfn(p, w, e)
xs = (1:length(w))';
expvals = p(1) + p(2) * exp(-xs / p(3));
err = sum((w - expvals).^2 ./ e);

if nargout >= 2
    gradA = -2 * sum((w - expvals) ./ e);
    gradB = 2 * sum((w - expvals) .* exp(-xs / p(3)) ./ e);
    gradTau = 2 * sum((w - expvals) ./ e .* (p(2) * exp(-xs / p(3)) / p(3)^2));
    grad = [gradA gradB gradTau];
end
end