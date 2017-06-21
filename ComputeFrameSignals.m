function [frame_signals] = ComputeFrameSignals(GaborData, oriKappa, spFreqCPP, spFreqStdCPP)
%COMPUTEFRAMESIGNALS wrapper around bpg.getSignal that computes the full
%[trials x frames] matrix of signal levels for the given data struct.
%
% If oriKappa is 0, the per-trial 'true' value of kappa is used. 
%
% spFreqCPP and spFreqStdCPP arguments are optional.

trials = length(GaborData.choice);

frame_signals = zeros(trials, GaborData.number_of_images);

if nargin == 2
    args = {oriKappa};
elseif nargin == 4
    args = {oriKappa, spFreqCPP, spFreqStdCPP};
else
    error('Expected 2 or 4 args');
end

parfor t=1:trials
    stim_images = GaborStimulus(GaborData, t);
    args_copy = args;
    if oriKappa == 0
        args_copy{1} = max(0.04, GaborData.noise(t));
    end
    frame_signals(t, :) = ...
        bpg.getSignal(stim_images, GaborData.left_category, args_copy{:}) - ...
        bpg.getSignal(stim_images, GaborData.right_category, args_copy{:});
end

end