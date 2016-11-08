function H = estEntropy(p, range, n)
%ESTENTROPY(p, [xmin xmax], [n]) estimate the entropy in the function
%handle p in the given range using bins defined by linspace(xmin, xmax, n).
%n is 10000 by default.
if nargin < 3, n = 10000; end
xs = linspace(range(1), range(2), n);
binwidth = xs(2)-xs(1);
pvals = arrayfun(p, xs);
pvals = pvals(pvals > 0);
H = -dot(pvals * binwidth, log(pvals));
end