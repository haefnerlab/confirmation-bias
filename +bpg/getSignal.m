function sig = getSignal(bpg_im, oriDEG, oriKappa, spFreqCPP, spFreqStdCPP)
%BPG.GETSIGNAL(bpg_im, oriDEG, oriKappa) compute the amount of energy there
%is in the given image(s) at the given orientation band. Spatial Frequency
%is integrated out.
%
%BPG.GETSIGNAL(bpg_im, oriDEG, oriKappa, spFreqCPP, spFreqStdCPP) also
%includes a spatial frequency filter.
%
% bpg_im must be an image (a matrix) or a [frames x height x width] array
% of images.
%
% Returns a [frames x 1] vector of signal levels.

if ismatrix(bpg_im)
    bpg_im = reshape(bpg_im, [1 size(bpg_im)]);
end
[frames, sz, ~] = size(bpg_im);

%% Construct orientation-band filter for the fourier domain.

[~, ~, rho, theta] = freq_coords(sz);
fouriFilter = bpg.vmpdf(2 * theta, 2 * deg2rad(oriDEG), oriKappa);
fouriFilter(isnan(fouriFilter)) = 0;

%% Include spatial frequency filter if specified

if nargin >= 5
    spFreqFilter = pdf('rician', rho / sz, spFreqCPP, spFreqStdCPP);
    fouriFilter = fouriFilter .* spFreqFilter;
end

%% Normalize the filter

fouriFilter = fouriFilter / sum(fouriFilter(:));

%% Get signal of each frame.
sig = zeros(frames, 1);
for f=1:frames
    frameF = fftshift(fft2(squeeze(bpg_im(f, :, :))));
    sig(f) = dot(fouriFilter(:), abs(frameF(:)));
end
end