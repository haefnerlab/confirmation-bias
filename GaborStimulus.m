function [image_array, frame_categories, checksum] = GaborStimulus(GaborData, trial)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial. If 'GaborData.iid(trial)' is true, each frame's category is
%drawn iid based on the 'ratio' parameter. Otherwise, exactly
%round(ratio*num_images) frames will match the 'true' category.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
rng(GaborData.seed(trial), 'twister');

if isfield(GaborData, 'flag_use_old_stimulus_code') && GaborData.flag_use_old_stimulus_code
    stim_fcn = @bpg.genImagesOld;
else
    stim_fcn = @bpg.genImages;
end

if ~isfield(GaborData, 'iid') || GaborData.iid(trial)
    % Randomly set each frame to match (or mismatch) the correct choice
    % for this trail, using the current 'ratio' to decide.
    match_frames = rand(1, GaborData.number_of_images) <= GaborData.ratio(trial);
else
    % Randomly permute whether each frame matches the true category, with
    % 'ratio' percent of them matching.
    n_match = round(GaborData.ratio(trial) * GaborData.number_of_images);
    match_frames = [true(1, n_match) false(1, GaborData.number_of_images - n_match)];
    match_frames = logical(Shuffle(double(match_frames)));
end

frame_categories = zeros(size(match_frames));

% Choose frames based on whether correct answer this trial is Left or Right
if GaborData.correct_answer(trial) == 1
    frame_categories(match_frames) = GaborData.left_category;
    frame_categories(~match_frames) = GaborData.right_category;
else
    frame_categories(~match_frames) = GaborData.left_category;
    frame_categories(match_frames) = GaborData.right_category;
end

% Set random seed again to keep match_frames independent of pixel noise.
rng(GaborData.seed(trial), 'twister');
image_array = stim_fcn(GaborData.number_of_images, GaborData.stim_size, ...
    GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, ...
    frame_categories, GaborData.noise(trial), GaborData.annulus);

image_array = uint8(image_array * GaborData.contrast(trial) + 127);

image_array = min(image_array, 255);
image_array = max(image_array, 0);

checksum = mean(image_array(:));

if ~GaborData.iid(trial) && (~isfield(GaborData, 'checksum') || GaborData.checksum(trial) == 0)
    warning(['Due to oversight, re-generating stimuli from seed *likely* fails on trials where iid=0. ' ...
        'Without a ''checksum'' field, there is no recovering!!']);
end

if isfield(GaborData, 'checksum') && GaborData.checksum(trial) ~= 0 && GaborData.checksum(trial) ~= checksum
    % Before erroring, try a hack-y solution... IF this trial had 'iid' set to 0, then it called the
    % Shuffle() function above, which has internal pseudorandomness completely separate from the RNG
    % state. Chances are, the current failure to recreate the stimulus is because Shuffle() changed,
    % not because the seed was wrong. Here, we catch that case by trying all permutations of the
    % shuffle.
    if ~GaborData.iid(trial)
        all_shuffles = unique(perms(match_frames), 'rows');
        for i=1:size(all_shuffles, 1)
            match_frames = all_shuffles(i, :);
            if GaborData.correct_answer(trial) == 1
                frame_categories(match_frames) = GaborData.left_category;
                frame_categories(~match_frames) = GaborData.right_category;
            else
                frame_categories(~match_frames) = GaborData.left_category;
                frame_categories(match_frames) = GaborData.right_category;
            end
            
            % Set random seed again to keep match_frames independent of pixel noise.
            rng(GaborData.seed(trial), 'twister');
            image_array = stim_fcn(GaborData.number_of_images, GaborData.stim_size, ...
                GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, ...
                frame_categories, GaborData.noise(trial), GaborData.annulus);
            
            image_array = uint8(image_array * GaborData.contrast(trial) + 127);
            
            image_array = min(image_array, 255);
            image_array = max(image_array, 0);
            
            checksum = mean(image_array(:));
            if checksum == GaborData.checksum(trial)
                warning('Checksum failed, but recovered by searching frame permutations...');
                return;
            end
        end
    end
    error('Stimulus reconstruction checksum failed!');
end

end