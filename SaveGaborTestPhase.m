function [] = SaveGaborTestPhase(subjectID, phase, loops, preliminary_trials, test_phase, directory)

% Save a subsection of the data from the experiment as test data manually
% chosen instead of automatically doing so in the experiment

if phase == 0       % load the test contrast data
    filename = fullfile(directory, 'RawData', [subjectID '-GaborDataContrast.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
elseif phase == 1       % load the test ratio data
    filename = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
end


test_start = 1;
test_end = 0;
Test_Data.current_trial = 0;
for t = 1:loops
    stage = test_phase(t);
    test_end =  test_end + preliminary_trials*t - stage + 1;
    
    Test_Data.move_on(test_start:test_end) = Preliminary_Data.move_on(stage:preliminary_trials*t);
    Test_Data.step_size(test_start:test_end) = Preliminary_Data.step_size(stage:preliminary_trials*t);
    Test_Data.reversal_counter(test_start:test_end) = Preliminary_Data.reversal_counter(stage:preliminary_trials*t);
    Test_Data.contrast(test_start:test_end) = Preliminary_Data.contrast(stage:preliminary_trials*t);
    Test_Data.correct_answer(test_start:test_end) = Preliminary_Data.correct_answer(stage:preliminary_trials*t);
    Test_Data.reaction_time(test_start:test_end) = Preliminary_Data.reaction_time(stage:preliminary_trials*t);
    Test_Data.choice(test_start:test_end) = Preliminary_Data.choice(stage:preliminary_trials*t);
    Test_Data.accuracy(test_start:test_end) = Preliminary_Data.accuracy(stage:preliminary_trials*t);
    Test_Data.ratio(test_start:test_end) = Preliminary_Data.ratio(stage:preliminary_trials*t);
    Test_Data.current_trial = Test_Data.current_trial + (preliminary_trials*t - stage + 1);
    Test_Data.test_phase = test_phase;
    Test_Data.number_of_images = Preliminary_Data.number_of_images;
    Test_Data.log_odds(test_start:test_end) = Preliminary_Data.log_odds(stage:preliminary_trials*t);
    Test_Data.average_orientations(:,test_start:test_end) = Preliminary_Data.average_orientations(:,stage:preliminary_trials*t);
    Test_Data.order_of_orientations(test_start:test_end,:) = Preliminary_Data.order_of_orientations(stage:preliminary_trials*t,:);
    Test_Data.screen_frame = Preliminary_Data.screen_frame;
    Test_Data.screen_resolution = Preliminary_Data.screen_resolution;
    Test_Data.image_length_x = Preliminary_Data.image_length_x;
    Test_Data.image_length_y = Preliminary_Data.image_length_y;
    Test_Data.image_template1(test_start:test_end,:) = Preliminary_Data.image_template1(stage:preliminary_trials*t,:);
    Test_Data.image_template2(test_start:test_end,:) = Preliminary_Data.image_template2(stage:preliminary_trials*t,:);
    Test_Data.image_template_difference(test_start:test_end,:) = Preliminary_Data.image_template_difference(stage:preliminary_trials*t,:);
    
    test_image_collection(test_start:test_end,:,:,:) = image_collection(stage:preliminary_trials*t,:,:,:);
    
    test_start = test_end+1;
end



if phase == 0       % save the test contrast data
    %% Save final data to folder
    if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
        mkdir(fullfile(directory, 'RawData'));
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestContrastManual.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data', 'test_image_collection'); % save the data
    else
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestContrastManual.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data', 'test_image_collection'); % save the data
    end
    
elseif phase == 1       % save the test ratio data
    %% Save final data to folder
    if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
        mkdir(fullfile(directory, 'RawData'));
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestRatioManual.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data', 'test_image_collection'); % save the data
    else
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestRatioManual.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data', 'test_image_collection'); % save the data
    end
end

end