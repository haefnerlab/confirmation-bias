function [pmf, lo, hi] = propagateNoisePMF(x_edges, pmf, shift, scale, noise)
assert(length(x_edges) == size(pmf,2)+1, 'x_edges expected to have size 1 more than elements of pmf');
assert(issorted(x_edges));

%% Stage 1: linearly remap x
% Begin by converting probability mass (pmf) into cumulative distribution (cmf) whose points
% correspond to x_edges
cmf = [0 cumsum(pmf, 2)];
remap_x = x_edges .* scale + shift;

if remap_x(end) < x_edges(1) || interp1(remap_x, cmf, x_edges(1)) > 0.999999
    % Edge case 1: everything is off the left edge. WARNING: known bug here: doesn't account for
    % possibility of all mass going out-of-bounds then coming back due to noise in stage 2.
    lo = 1;
    hi = 0;
    pmf = zeros(size(pmf));
    return;
elseif remap_x(1) > x_edges(end) || interp1(remap_x, cmf, x_edges(end)) < 0.000001
    % Edge case 2: everything is off the right edge. See warning in previous block.
    lo = 0;
    hi = 1;
    pmf = zeros(size(pmf));
    return;
end

% Having cleared the above two edge cases, we can guarantee at least one non-nan value in remap_cmf.

remap_cmf = interp1(remap_x, cmf, x_edges, 'linear');
firstidx = find(~isnan(remap_cmf), 1);
if firstidx, remap_cmf(1:firstidx-1) = 0; end
lastidx = find(~isnan(remap_cmf), 1, 'last');
if lastidx, remap_cmf(lastidx+1:end) = 1; end

% Count total mass that went out-of bounds (lo or hi)
lo = interp1(remap_x, cmf, x_edges(1), 'linear');
if isnan(lo) || lo < 0, lo = 0; end
hi = interp1(remap_x, 1-cmf, x_edges(end), 'linear');
if isnan(hi) || hi < 0, hi = 0; end

% Ensure valid CMF in the "middle" section
remap_cmf = remap_cmf - remap_cmf(1);
remap_cmf = remap_cmf / remap_cmf(end);

pmf = diff(remap_cmf) * (1 - (lo + hi));

if any(isnan(remap_cmf))
    keyboard
end

if noise == 0, return; end

%% Stage 2: add noise
dx = diff(x_edges);
mdx = mean(dx);
if all(abs(dx-mdx) < min(dx)*1e-6)
    % Uniform grid: adding noise is convolution
    bins_per_5sig = ceil(5*noise/mdx);
    noise_x_edges = linspace(-mdx*bins_per_5sig, mdx*bins_per_5sig, 2*(bins_per_5sig+1));
    noise_x_center = (noise_x_edges(1:end-1)+noise_x_edges(2:end))/2;
    noise_pmf = normpdf(noise_x_center, 0, noise);
    noise_pmf = noise_pmf / sum(noise_pmf);
    noise_pmf = (noise_pmf + fliplr(noise_pmf)) / 2;
    
    pmf_bleed = conv(pmf, noise_pmf, 'full');
    lo = lo + sum(pmf_bleed(1:bins_per_5sig));
    hi = hi + sum(pmf_bleed(end-bins_per_5sig+1:end));
    pmf = pmf_bleed(bins_per_5sig+1:end-bins_per_5sig);
else
    % Non-uniform grid: loop over each x to add noise
    error('non-unform grid not yet implemented');
end

total_mass = lo + hi + sum(pmf);
assert(abs(total_mass - 1) < 1e-9, 'wtf');

end