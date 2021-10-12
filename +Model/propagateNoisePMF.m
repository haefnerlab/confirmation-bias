function [pmf, lo, hi] = propagateNoisePMF(x_edges, pmf, shift, scale, noise)
assert(size(x_edges,2) == size(pmf,2)+1, 'x_edges expected to have size 1 more than elements of pmf');
assert(issorted(x_edges,2));

N = size(pmf,1);
if size(x_edges,1) == 1 && N > 1
    x_edges = repmat(x_edges, N, 1);
end

%% Stage 1: linearly remap x
% Begin by converting probability mass (pmf) into cumulative distribution (cmf) whose points
% correspond to x_edges
cmf = [zeros(N,1) cumsum(pmf, 2)];
remap_x = x_edges .* scale + shift;

% Find rows where *all* data went out of the bounds of x_edges. These will fail interpolation in the
% next step.
all_oob_lo = remap_x(:,end) < x_edges(:,1) | myinterp1(remap_x, cmf, x_edges(:,1)) > 0.999999;
all_oob_hi = remap_x(:,1) > x_edges(:,end) | myinterp1(remap_x, cmf, x_edges(:,end)) < 0.000001;

% Get new CMF in old x_edges coordinates by linear interpolation of the shifted CMF. By default,
% indices of x_edges that are out of bounds should get the correct '0' or '1' CMF value, unless ALL
% indices are out of bounds, but that is handled by 'all_oob_lo' and 'all_oob_hi'
remap_cmf = myinterp1(remap_x, cmf, x_edges);

% Count total mass that went out-of bounds (lo or hi)
lo = myinterp1(remap_x, cmf, x_edges(:,1));
lo(isnan(lo) | lo < 0) = 0;
hi = myinterp1(remap_x, 1-cmf, x_edges(:,end));
hi(isnan(hi) | hi < 0) = 0;

% Ensure valid CMF in the "middle" section
remap_cmf = remap_cmf - remap_cmf(:,1);
remap_cmf = remap_cmf ./ remap_cmf(:,end);

pmf = diff(remap_cmf,[],2) .* (1 - (lo + hi));

% Handle the 'all out of bounds' cases
pmf(all_oob_lo,:) = 0;
pmf(all_oob_hi,:) = 0;
lo(all_oob_lo) = 1;
hi(all_oob_lo) = 0;
lo(all_oob_hi) = 0;
hi(all_oob_hi) = 1;

if noise == 0, return; end

%% Stage 2: add noise.
dx = diff(x_edges,[],2);
mdx = mean(dx(:));
assert(all(abs(dx(:) - mdx) < 1e-9), 'Cannot handle unqually spaced grids!');

% Uniform grid: adding noise is convolution. Make kernel out to 5 sigma of noise.
bins_per_5sig = ceil(5*noise/mdx);
noise_x_edges = linspace(-mdx*bins_per_5sig, mdx*bins_per_5sig, 2*(bins_per_5sig+1));
noise_x_center = (noise_x_edges(1:end-1)+noise_x_edges(2:end))/2;
noise_pmf = normpdf(noise_x_center, 0, noise);
noise_pmf = noise_pmf / sum(noise_pmf);
noise_pmf = (noise_pmf + fliplr(noise_pmf)) / 2;

pmf_bleed = conv2(1, noise_pmf, pmf, 'full');
lo = lo + sum(pmf_bleed(:, 1:bins_per_5sig),2);
hi = hi + sum(pmf_bleed(:, end-bins_per_5sig+1:end),2);
pmf = pmf_bleed(:, bins_per_5sig+1:end-bins_per_5sig);

total_mass = lo + hi + sum(pmf, 2);
assert(all(abs(total_mass - 1) < 1e-9), 'wtf');

end

function val = myinterp1(x, v, xq)
% Faster linear interpolation (same signature as interp1) assuming x is an evenly spaced and sorted
% row vector. May be run vectorized where x and v are [N x grid] and xq is [N x k] for k queries on
% each of N grids. Out-of-bounds queries inherit the value of v(1) or v(end).

[N, g] = size(x);
assert(all(size(x) == size(v)), 'Size of x and v must match');
assert(size(x,1) == size(xq,1), 'Number of rows of x and xq must match');

% Detect 'out of bounds' queries. These will be set to nan.
oob_left = xq <= x(:,1);
oob_right = xq  >= x(:,end);

% Find the index into 'grid' for each query such that x(idx) < xq < x(idx+1)
grid_frac = (xq - x(:,1)) ./ (x(:,end)-x(:,1));
idx = floor((g-1) * grid_frac) + 1;

% Avoid index errors for out-of-bounds xq values
idx(oob_left | oob_right) = 1;

% Convert from 2D into 1D index for vectorized case.
rowno = (1:N)';
idx1d_lo = (idx-1)*N+rowno;
idx1d_hi = idx*N+rowno;

% Interpolation weight ranges in [0 1] for each xq specifying the fraction of the way between x(idx)
% and x(idx+1) it lies.
wt = (xq-x(idx1d_lo))./(x(idx1d_hi)-x(idx1d_lo));
val = v(idx1d_hi).*wt + v(idx1d_lo).*(1-wt);

val(oob_left) = v(1);
val(oob_right) = v(end);
end