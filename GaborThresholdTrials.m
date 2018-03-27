function Data = GaborThresholdTrials(Data, phase, threshold, floor)

if phase == 0
    stair_param = 'contrast';
elseif phase == 1
    stair_param = 'true_ratio';
elseif phase == 2
    stair_param = 'noise';
end

test_trials = Data.(stair_param) <= threshold;
if exist('floor', 'var')
    test_trials = test_trials & Data.(stair_param) >= floor;
end
elements = sum(test_trials);

%% Subselect trial data.

Data.current_trial = elements;

Data.contrast = Data.contrast(test_trials);
Data.ratio = Data.ratio(test_trials);
Data.noise = Data.noise(test_trials);
Data.step_size = Data.step_size(test_trials);

Data.iid = Data.iid(test_trials);
Data.seed = Data.seed(test_trials);
if isfield(Data, 'checksum'), Data.checksum = Data.checksum(test_trials); end
Data.streak = Data.streak(test_trials);
Data.reversal_counter = Data.reversal_counter(test_trials);
Data.correct_answer = Data.correct_answer(test_trials);
Data.ideal_answer = Data.ideal_answer(test_trials);
Data.reaction_time = Data.reaction_time(test_trials);
Data.choice = Data.choice(test_trials);
Data.accuracy = Data.accuracy(test_trials);
Data.frame_categories = Data.frame_categories(test_trials, :);
Data.ideal_frame_signals = Data.ideal_frame_signals(test_trials, :);

if isempty(Data.model_observer)
    Data.eye_tracker_points = Data.eye_tracker_points(test_trials);
end

end
