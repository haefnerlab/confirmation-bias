function [images, signals] = regenSignals(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, contrast, annulusPix, filterKappa)
%% Get signal on generated images
[images, ~] = bpg.genImages(nSamples, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix);
images = uint8(images * contrast + 127);
images = min(images, 255);
images = max(images, 0);

images = double(images)-127;
signals = bpg.getSignal(images, oriDEG, filterKappa) - bpg.getSignal(images, oriDEG-90, filterKappa);
end