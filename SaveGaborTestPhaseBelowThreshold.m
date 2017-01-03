function [Test_Data, image_collection_test] = SaveGaborTestPhaseBelowThreshold(subjectID, phase, threshold, directory)

cd(fullfile(directory, 'Code'));
savedir = fullfile(directory, 'RawData');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if phase == 0
    postfix = '-GaborDataContrast.mat';
    postfix_test = '-GaborTestContrast.mat';
    stair_param = 'contrast';
elseif phase == 1
    postfix = '-GaborDataRatio.mat';
    postfix_test = '-GaborTestRatio.mat';
    stair_param = 'ratio';
end

filename = fullfile(directory, 'RawData', [subjectID postfix]);
if ~exist(filename, 'file')
    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
    disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
    return;
else
    contents = load(filename); % Load Preliminary_Data and image_collection
end

Data = contents.GaborData;
image_collection = contents.image_collection;

test_trials = Data.(stair_param) < threshold;
elements = sum(test_trials);

%% Copy from image_collection to image_collection_test

image_collection_test = image_collection(test_trials, :, :, :);

%% Copy from Data to Test_Data

Test_Data = Data;

Test_Data.current_trial = elements;

Test_Data.contrast = Data.contrast(test_trials);
Test_Data.ratio = Data.ratio(test_trials);
Test_Data.pixel_noise = Data.pixel_noise(test_trials);
Test_Data.step_size = Data.step_size(test_trials);

Test_Data.streak = Data.streak(test_trials);
Test_Data.reversal_counter = Data.reversal_counter(test_trials);
Test_Data.correct_answer = Data.correct_answer(test_trials);
Test_Data.ideal_answer = Data.ideal_answer(test_trials);
Test_Data.reaction_time = Data.reaction_time(test_trials);
Test_Data.choice = Data.choice(test_trials);
Test_Data.accuracy = Data.accuracy(test_trials);
Test_Data.order_of_orientations = Data.order_of_orientations(test_trials, :);
Test_Data.log_frame_odds = Data.log_frame_odds(test_trials, :);
Test_Data.log_decision_odds = Data.log_decision_odds(test_trials, :);
Test_Data.average_orientations = Data.average_orientations(:, test_trials);

Test_Data.eye_tracker_points = Data.eye_tracker_points(test_trials);

%% Save Test_Data and image_collection_test

savefile = fullfile(savedir, [subjectID postfix_test]);
save(savefile, 'Test_Data', 'image_collection_test');

end
