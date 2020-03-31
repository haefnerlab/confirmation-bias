function [image_array, frame_categories, checksum] = GaborStimulus(GaborData, trial, regen)
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
    idx_mismatch = randperm(GaborData.number_of_images, GaborData.number_of_images - n_match);
    match_frames = true(1, GaborData.number_of_images);
    match_frames(idx_mismatch) = false;
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

%% Done generating stim... perform tests to ensure accuracy if regenerating from seed
if nargin >= 3 && regen > 0
    % We are regenerating from seed... if regen == 0, just blindly return. If regen == 1, comare to
    % checksum field (errors are rare but possible). If regen == 2, compare to ideal_frame_signals
    % field (slower but less error prone).
    if regen == 1 && (~isfield(GaborData, 'checksum') || GaborData.checksum(trial) == 0)
        warning('GABORSTIMULUS : flag regen = 1, but no available checksum information! Defaulting to regen = 2');
        regen = 2;
    end
    
    if regen == 1
        % Do verification based on checksum
        if checksum == GaborData.checksum(trial)
            % Simple enough... checksum passed (BUT errors may nonetheless slip through)
            return
        else
            % Checksum failed. Before erroring, try a hack-y solution... IF this trial had 'iid' set
            % to 0, then it called the Shuffle() function above, which has internal pseudorandomness
            % completely separate from the RNG state. Chances are, the current failure to recreate
            % the stimulus is because Shuffle() changed, not because the seed was wrong. Here, we
            % catch that case by trying all permutations of the shuffle.
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
                        warning('Checksum-by-sum failed, but recovered by searching frame permutations...');
                        return;
                    end
                end
            end
            error('Stimulus reconstruction checksum failed (by checksum field)!');
        end
    elseif regen == 2
        % Do verification based on ideal_frame_signals. NOTE: image_array is type uint8, so it is an
        % integer in [0 255]. Subtracting 127 clips values below zero, i.e. the result is in [0
        % 128]. This is a bug, but it replicates an identical bug in @ExperimentGabor when
        % populating the 'ideal_frame_signals' field. Our goal here is to regenerate the identical
        % signals to verify that the stimulus was reconstructed properly, not to get the 'true'
        % signal values per se. (Also, testing has shown that despite the zero clipping here, the
        % 'true' signals (i.e. calling image_array=double(image_array) first) are very highly
        % correlated with these values.)
        stored_signals = GaborData.ideal_frame_signals(trial, :);
        kernel_kappa = max(GaborData.noise(trial), .04);
        regen_signals = bpg.getSignal(image_array-127, GaborData.left_category, kernel_kappa) - ...
            bpg.getSignal(image_array-127, GaborData.right_category, kernel_kappa);
        
        if all(abs(regen_signals' - stored_signals) < 1e-10)
            % Passed!
            return
        else
            % See above comment on the 'iid' flag...
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
                    
                    regen_signals = bpg.getSignal(image_array-127, GaborData.left_category, kernel_kappa) - ...
                        bpg.getSignal(image_array-127, GaborData.right_category, kernel_kappa);
                    if all(abs(regen_signals' - stored_signals) < 1e-10)
                        warning('Checksum-by-signal failed, but recovered by searching frame permutations...');
                        return;
                    end
                end
            end
            error('Stimulus reconstruction checksum failed (by signals)!');
        end
    end
end
end