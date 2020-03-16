function [images, signals] = regenSignals(nSamples, nFramesPerTrial, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, contrast, annulusPix, filterKappa)
%% Get signal on generated images

% Round-up # samples to nearest whole trial
nTrials = ceil(nSamples / nFramesPerTrial);
nSamples = nFramesPerTrial * nTrials;

[images, ~] = bpg.genImages(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix);

% Important: normalize to max in batches of nFramesPerTrial at a time. This is important because
% bpg.genImages normalizes by the MAX pixel value of everything generated. When generating 100s or
% 1000s of frames at once, the max is likely to be very outlying. To recreate the distribution of
% signals on the actual data, we need to normalize to the max of each batch of 10.
max_per_image = max(abs(reshape(images, nSamples, width*width)), [], 2);
% Group images into 'nFramesPerTrial' and compute the max pixel value per pseudo-trial
max_per_trial = max(reshape(max_per_image, nTrials, nFramesPerTrial), [], 2);
% Normalize imagesby max pixel value trial-wise 
max_norm = repmat(max_per_trial, 1, nFramesPerTrial);
images = images ./ max_norm(:);

% Copy of uint8 conversion that appears in @GaborStimulus
images = uint8(images * contrast + 127);
images = min(images, 255);
images = max(images, 0);

% Copy of centering and signal computation that appears in @ComputeFrameSignals
images = double(images)-127;
signals = bpg.getSignal(images, oriDEG, filterKappa) - bpg.getSignal(images, oriDEG-90, filterKappa);
end