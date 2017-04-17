function [im, imF] = genImages(frames, width, spFreqCPP, spFreqStdCPP, oriDEG, oriStdDEG)
%BPG.GENIMAGES Create a sequence band-pass grating (bpg) stimuli.
%
%[im, imF] = BPG.GENIMAGES(frames, width, spFreqCPP, spFreqStdCPP, oriDEG, oriStdDEG)
% creates [frames x width x width] array of images. spFreqCPP sets the mean
% spatial frequency in cycles per pixel. spFreqStdCPP sets the range of
% spatial frequencies present. oriDEG sets the mean rotation, oriStdDeg
% sets the range of orientation energy present.

noise = randn(frames, width, width);
noiseF = framefun(@(f) fftshift(fft2(f)), noise);

x = linspace(-1, 1, width);
[xx, yy] = meshgrid(x, x);
% Get polar coordinates: rr is radius, tt is angle.
rr = sqrt(xx.^2 + yy.^2);
tt = atan2(yy, xx);

%% Create spatial frequency filter
spFreqFilter = pdf('rician', rr, spFreqCPP, spFreqStdCPP);

%% Create orientation filter
oriFilter = bpg.vmpdf(2 * tt, deg2rad(oriDEG), 1 / deg2rad(oriStdDEG));

%% Apply fourier-domain filters
filterF = spFreqFilter .* oriFilter;
filterF = filterF / max(filterF(:));

imF = framefun(@(f) f .* filterF, noiseF);
im = framefun(@(f) real(ifft2(ifftshift(f))), imF);

%% Create gaussian aperture
aperture = exp(-4 * rr.^2);

im = framefun(@(f) f .* aperture, im);

end

function frames = framefun(fn ,frames)
%Helper to apply fn to each frame in frames.
for f=1:size(frames, 1)
    frames(f, :, :) = fn(squeeze(frames(f, :, :)));
end
end