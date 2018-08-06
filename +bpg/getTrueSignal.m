function [im, imF, sig, filterF, aperture] = getTrueSignal(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix, filterKappa)

noise = randn(width, width, nSamples);
for iSample=nSamples:-1:1
    noiseF(:,:,iSample) = fftshift(fft2(noise(:,:,iSample)));
    noiseF(:,:,iSample) = noiseF(:,:,iSample) ./ abs(noiseF(:,:,iSample));
end

[~, ~, rho, theta] = freq_coords(width);

% Create separate [-1, 1] range meshgrid for pixel-space filters.
[px, py] = meshgrid(linspace(-1, 1, width));
pr = sqrt(px.^2 + py.^2);

%% Create spatial frequency filter
spFreqFilter = pdf('rician', rho / width, spFreqCPP, spFreqStdCPP);

%% Create gaussian aperture
aperture = exp(-4 * pr.^2);

if nargin >= 7 && annulusPix > 0
    % Cut out annulus hole.
    aperture = aperture .* (1 + erf(10 * (pr - annulusPix / width)));
end

%% Apply filters and get image

% Create orientation filters for each frame. Note that 'theta' is
% doubled to create two symmetric filters in the Fourier domain
% (bow-tie rather than cone shape). 'oriDEG' must also be doubled to
% compensate.
oriFilter = bpg.vmpdf(2 * theta, 2 * deg2rad(oriDEG), oriKappa);
oriFilter(isnan(oriFilter)) = 0;

% Get full, normalized foureir-domain filter.
filterF = spFreqFilter .* oriFilter;
filterF = filterF / sum(filterF(:));

% Apply fourier-domain filters.
imF = noiseF .* filterF;
for iSample=nSamples:-1:1
    im(:,:,iSample) = aperture .* real(ifft2(ifftshift(imF(:,:,iSample))));
end

%% Normalize range in pixel space to +/- 1
im = im / max(abs(im(:)));

%% Get signal on generated image

im = permute(im, [3 1 2]);
sig = bpg.getSignal(im, oriDEG, filterKappa) - bpg.getSignal(im, oriDEG+90, filterKappa);
end