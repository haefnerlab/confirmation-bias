function GaborData = ratio(GaborData)
%RATIO updates staircase-able parameters (contrast, ratio, and
%pixel_noise) in GaborData for the current trial based on the results of
%the previous trial. Only ratio may change to alter the difficulty.
%Requires the following are already set appropriately:
%
% GaborData.contrast(trial-1)
% GaborData.ratio(trial-1)
% GaborData.pixel_noise(trial-1)
% GaborData.step_size(trial-1)
% GaborData.streak(trial)
% GaborData.reversal_counter(trial)
%
%Computes and sets the following:
%
% GaborData.contrast(trial)
% GaborData.ratio(trial)
% GaborData.pixel_noise(trial)
% GaborData.step_size(trial)

trial = GaborData.current_trial;

%% Copy over params - contrast may be overwritten below
GaborData.contrast(trial) = GaborData.contrast(trial-1);
GaborData.ratio(trial) = GaborData.ratio(trial-1);
GaborData.pixel_noise(trial) = GaborData.pixel_noise(trial-1);

%% Reduce step size after 10 reversals
prev_reversals = GaborData.reversal_counter(trial-1);
reversals = GaborData.reversal_counter(trial);
m = GaborData.reversals_per_epoch;
if reversals > 0 && mod(reversals, m) == 0 && mod(prev_reversals, m) ~= 0
    % Decay the step size half way towards 0
    GaborData.step_size(trial) = GaborData.step_size(trial-1) / 2;
else
    % Same step size as last trial
    GaborData.step_size(trial) = GaborData.step_size(trial-1);
end

if GaborData.step_size(trial) < GaborData.min_step_size
    GaborData.step_size(trial) = GaborData.min_step_size;
end

%% Apply staircase logic
if GaborData.streak(trial) == 0
    % Got it wrong - make things easier
    GaborData.ratio(trial) = ...
        GaborData.ratio(trial-1) + GaborData.step_size(trial);
elseif mod(GaborData.streak(trial), 2) == 0
    % Got 2 right in a row - make things harder
    GaborData.ratio(trial) = ...
        GaborData.ratio(trial-1) - GaborData.step_size(trial);
end

% Apply bounds
GaborData.ratio(trial) = ...
    max(GaborData.ratio(trial), GaborData.stair_bounds(1));
GaborData.ratio(trial) = ...
    min(GaborData.ratio(trial), GaborData.stair_bounds(2));

end