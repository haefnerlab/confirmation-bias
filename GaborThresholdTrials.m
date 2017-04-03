function [Data, image_collection] = GaborThresholdTrials(Data, image_collection, phase, threshold, floor)

if phase == 0
    stair_param = 'contrast';
elseif phase == 1
    stair_param = 'ratio';
end

test_trials = Data.(stair_param) < threshold;
if exist('floor', 'var')
    test_trials = test_trials & Data.(stair_param) > floor;
end
elements = sum(test_trials);

%% Subselect images.

image_collection = image_collection(test_trials, :, :, :);

%% Subselect trial data.

Data.current_trial = elements;

Data.contrast = Data.contrast(test_trials);
Data.ratio = Data.ratio(test_trials);
Data.pixel_noise = Data.pixel_noise(test_trials);
Data.step_size = Data.step_size(test_trials);

Data.streak = Data.streak(test_trials);
Data.reversal_counter = Data.reversal_counter(test_trials);
Data.correct_answer = Data.correct_answer(test_trials);
Data.ideal_answer = Data.ideal_answer(test_trials);
Data.reaction_time = Data.reaction_time(test_trials);
Data.choice = Data.choice(test_trials);
Data.accuracy = Data.accuracy(test_trials);
Data.order_of_orientations = Data.order_of_orientations(test_trials, :);
Data.log_frame_odds = Data.log_frame_odds(test_trials, :);
Data.log_decision_odds = Data.log_decision_odds(test_trials, :);
Data.average_orientations = Data.average_orientations(:, test_trials);

Data.eye_tracker_points = Data.eye_tracker_points(test_trials);

end
