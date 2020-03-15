function [im, imF, sig, filterF] = regenSignals(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, contrast, annulusPix, filterKappa)
%% Get signal on generated images
[im, ~, filterF] = bpg.genImages(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix);
im = uint8(im * contrast + 127);
imF = fftshift(fft2(im-127));

% NOTE: im is type uint8, so it is an integer in [0 255]. Subtracting 127 clips values below zero,
% i.e. the result is in [0 128]. This is a bug, but it replicates an identical bug in
% @ExperimentGabor. Our goal here is to regenerate a set of images 'im' and signals 'sig' that have
% identical statistics as values generated during the experiment. Offline testing has shown that
% despite the zero clipping here, the 'true' signals (i.e. calling im=double(im) first) are very
% highly correlated with these values.
sig = bpg.getSignal(im - 127, oriDEG, filterKappa) - bpg.getSignal(im - 127, oriDEG+90, filterKappa);
end