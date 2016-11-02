function [] = SaveAuditoryTestPhaseBelowThreshold(subjectID, phase, loops, preliminary_trials,threshold, directory)

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
elements=0;
if phase ==0
    for t = 1:loops*preliminary_trials
        if (Preliminary_Data.volume(t)<=threshold)
            elements=elements+1;
        end
    end
elseif phase ==1
    for t = 1:loops*preliminary_trials
        if (Preliminary_Data.ratio(t)<=threshold)
            elements=elements+1;
        end
    end
end

Test_Data.current_trial = 0;
Test_Data.move_on = zeros(1,elements);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
Test_Data.step_size = zeros(1,elements);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
Test_Data.reversal_counter = zeros(1,elements);   % How many trials has the subject got wrong? When to change the step size?
Test_Data.volume = zeros(1,elements);         % How loud the sound is, or the signal level
Test_Data.correct_answer = zeros(1,elements);         % What was the right ear/answer?
Test_Data.reaction_time = zeros(1,elements);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Test_Data.choice = zeros(1,elements);                 % Which ear did the subject choose? 1 = left, 0 = right
Test_Data.accuracy = zeros(1,elements);               % Did the subject make the right selection? 1 = yes, 0 = no
Test_Data.ratio = zeros(1,elements);  % Re ords underlying click rate
Test_Data.staircase_answer = zeros(1,elements);%stores answers based on actual click rates for staircase


% Set properties of a click
Test_Data.sampling_rate = 6000;   % Indirectly controls frequency
Test_Data.stimulus_duration = 2;           % How many seconds is each trial?
Test_Data.bins = Preliminary_Data.sampling_rate * Preliminary_Data.stimulus_duration;         % There needs to be an number of bins for the sound vector

Test_Data.click_rate = zeros(2,elements);          % How many clicks did each ear actually hear?
Test_Data.average_clicks = zeros(2,elements);      % Average number of clicks to be the mean of the poisson distribution
Test_Data.number_of_frames = 60 * Test_Data.stimulus_duration; % How many 'frames' or chances for a click per second in a trial?

Test_Data.order_of_clicks = zeros(elements, 2, Test_Data.number_of_frames);
if phase==0
    k=1;
    for i = 1:loops*preliminary_trials
        if (Preliminary_Data.volume(i)<=threshold)
            
            Test_Data.move_on(k) = Preliminary_Data.move_on(i);
            Test_Data.step_size(k) = Preliminary_Data.step_size(i);
            Test_Data.reversal_counter(k) = Preliminary_Data.reversal_counter(i);
            Test_Data.volume(k) = Preliminary_Data.volume(i);
            Test_Data.correct_answer(k) = Preliminary_Data.correct_answer(i);
            Test_Data.staircase_answer(k) = Preliminary_Data.staircase_answer(i);
            Test_Data.reaction_time(k) = Preliminary_Data.reaction_time(i);
            Test_Data.choice(k) = Preliminary_Data.choice(i);
            Test_Data.accuracy(k) = Preliminary_Data.accuracy(i);
            Test_Data.ratio(k) = Preliminary_Data.ratio(i);
            Test_Data.current_trial = k;
            Test_Data.sampling_rate = Preliminary_Data.sampling_rate;
            Test_Data.stimulus_duration = Preliminary_Data.stimulus_duration;
            Test_Data.bins = Preliminary_Data.bins;
            Test_Data.click_rate(:,k) = Preliminary_Data.click_rate(:,i);
            Test_Data.average_clicks(:,k) = Preliminary_Data.average_clicks(:,i);
            Test_Data.number_of_frames = Preliminary_Data.number_of_frames;
            Test_Data.order_of_clicks(k,:,:) = Preliminary_Data.order_of_clicks(i,:,:);
            
            k=k+1;
        end
    end
elseif phase==1
    k=1;
    for i = 1:loops*preliminary_trials
        if (Preliminary_Data.ratio(i)<=threshold)
            
            Test_Data.move_on(k) = Preliminary_Data.move_on(i);
            Test_Data.step_size(k) = Preliminary_Data.step_size(i);
            Test_Data.reversal_counter(k) = Preliminary_Data.reversal_counter(i);
            Test_Data.volume(k) = Preliminary_Data.volume(i);
            Test_Data.correct_answer(k) = Preliminary_Data.correct_answer(i);
            Test_Data.staircase_answer(k) = Preliminary_Data.staircase_answer(i);
            Test_Data.reaction_time(k) = Preliminary_Data.reaction_time(i);
            Test_Data.choice(k) = Preliminary_Data.choice(i);
            Test_Data.accuracy(k) = Preliminary_Data.accuracy(i);
            Test_Data.ratio(k) = Preliminary_Data.ratio(i);
            Test_Data.current_trial = k;
            Test_Data.sampling_rate = Preliminary_Data.sampling_rate;
            Test_Data.stimulus_duration = Preliminary_Data.stimulus_duration;
            Test_Data.bins = Preliminary_Data.bins;
            Test_Data.click_rate(:,k) = Preliminary_Data.click_rate(:,i);
            Test_Data.average_clicks(:,k) = Preliminary_Data.average_clicks(:,i);
            Test_Data.number_of_frames = Preliminary_Data.number_of_frames;
            Test_Data.order_of_clicks(k,:,:) = Preliminary_Data.order_of_clicks(i,:,:);
            
            k=k+1;
        end
    end
end

if phase == 0       % save the test contrast data
    %% Save final data to folder
    if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-AuditoryTestVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-AuditoryTestVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
elseif phase == 1       % save the test ratio data
    %% Save final data to folder
    if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-AuditoryTestRatio.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-AuditoryTestRatio.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
end

end