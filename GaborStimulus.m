function [image_array, frame_categories, true_category] = GaborStimulus(GaborData, trial)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
rng(GaborData.seed(trial), 'twister');

% Randomly set each frame to match (or mismatch) the correct choice
% for this trail, using the current 'ratio' to decide.
match_frames = 1 * (rand(1, GaborData.number_of_images) <= GaborData.ratio(trial));

% Choose whether correct answer this trial will be Left or Right
if rand < 0.5
    frame_categories = match_frames;
    true_category = 1;
else
    frame_categories = 1 - match_frames;
    true_category = 0;
end

% Set random seed again to keep match_frames independent of pixel noise.
rng(GaborData.seed(t), 'twister');
image_array = makeImages(GaborData.image_length_x, frame_categories, ...
    GaborData.left_template, GaborData.right_template, ...
    GaborData.contrast(trial), GaborData.noise(trial));

end