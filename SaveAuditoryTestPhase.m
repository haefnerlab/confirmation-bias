function [] = SaveAuditoryTestPhase(subjectID, phase, loops, preliminary_trials, test_phase, directory)

% Save a subsection of the data from the experiment as test data manually
% chosen instead of automatically doing so in the experiment

if phase == 0       % load the test volume data
    filename = [directory 'RawData/' subjectID '-AuditoryDataVolume.mat'];
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
elseif phase == 1       % load the test ratio data
    filename = [directory 'RawData/' subjectID '-AuditoryDataRatio.mat'];
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
    Test_Data.volume(test_start:test_end) = Preliminary_Data.volume(stage:preliminary_trials*t);
    Test_Data.correct_answer(test_start:test_end) = Preliminary_Data.correct_answer(stage:preliminary_trials*t);
    Test_Data.staircase_answer(test_start:test_end) = Preliminary_Data.staircase_answer(stage:preliminary_trials*t);
    Test_Data.reaction_time(test_start:test_end) = Preliminary_Data.reaction_time(stage:preliminary_trials*t);
    Test_Data.choice(test_start:test_end) = Preliminary_Data.choice(stage:preliminary_trials*t);
    Test_Data.accuracy(test_start:test_end) = Preliminary_Data.accuracy(stage:preliminary_trials*t);
    Test_Data.ratio(test_start:test_end) = Preliminary_Data.ratio(stage:preliminary_trials*t);
    Test_Data.current_trial = Test_Data.current_trial + (preliminary_trials*t - stage + 1);
    Test_Data.test_phase = test_phase;
    Test_Data.sampling_rate = Preliminary_Data.sampling_rate;
    Test_Data.stimulus_duration = Preliminary_Data.stimulus_duration;
    Test_Data.bins = Preliminary_Data.bins;
    Test_Data.click_rate(:,test_start:test_end) = Preliminary_Data.click_rate(:,stage:preliminary_trials*t);
    Test_Data.average_clicks(:,test_start:test_end) = Preliminary_Data.average_clicks(:,stage:preliminary_trials*t);
    Test_Data.number_of_frames = Preliminary_Data.number_of_frames;
    Test_Data.order_of_clicks(test_start:test_end,:,:) = Preliminary_Data.order_of_clicks(stage:preliminary_trials*t,:,:);
    
    test_start = test_end+1;
end



if phase == 0       % save the test contrast data
    %% Save final data to folder
    if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-AuditoryTestVolumeManual.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-AuditoryTestVolumeManual.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
elseif phase == 1       % save the test ratio data
    %% Save final data to folder
    if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-AuditoryTestRatioManual.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-AuditoryTestRatioManual.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
end

end