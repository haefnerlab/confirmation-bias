function expfit = expFit(weights, errors)
if nargin < 2, errors = ones(size(weights)); end

xs = (1:length(weights))';

expfn = @(x, p) p(1) + p(2) * exp(-x / p(3));

errfn = @(p) sum((weights - expfn(xs, p)).^2 ./ errors);

p0 = [mean(weights) 0 10];

expfit = fmincon(errfn, p0, [], [], [], [], [0 -inf 0], [inf inf inf]);

end