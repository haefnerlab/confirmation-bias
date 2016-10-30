function [] = SaveQuitDataAuditory(subjectID, phase, loops_consider, preliminary_trials, directory)

if phase == 0       % load the volume data
    filename = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
elseif phase == 1       % load the ratio data
    filename = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
end

point = loops_consider*preliminary_trials;
Preliminary_Data.move_on(point+1:end) = [];          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
Preliminary_Data.step_size(point+1:end) = [];        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
Preliminary_Data.reversal_counter(point+1:end) = [];   % How many trials has the subject got wrong? When to change the step size?
Preliminary_Data.volume(point+1:end) = [];         % How loud the sound is, or the signal level
Preliminary_Data.correct_answer(point+1:end) = [];         % What was the right ear/answer?
Preliminary_Data.reaction_time(point+1:end) = [];          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Preliminary_Data.choice(point+1:end) = [];                 % Which ear did the subject choose? 1 = left, 0 = right
Preliminary_Data.accuracy(point+1:end) = [];               % Did the subject make the right selection? 1 = yes, 0 = no
Preliminary_Data.ratio(point+1:end) = [];  % Re ords underlying click rate
Preliminary_Data.staircase_answer(point+1:end) = [];%stores answers based on actual click rates for staircase

Preliminary_Data.current_trial = point;
%Preliminary_Data.test_phase = ([1:loops_consider].*preliminary_trials) - preliminary_trials + 1;
Preliminary_Data.test_phase ((loops_consider+1):end) = [];

% Set properties of a click
Preliminary_Data.sampling_rate = 6000;   % Indirectly controls frequency
Preliminary_Data.stimulus_duration = 2;           % How many seconds is each trial?
Preliminary_Data.bins = Preliminary_Data.sampling_rate * Preliminary_Data.stimulus_duration;         % There needs to be an number of bins for the sound vector

Preliminary_Data.click_rate(:,point+1:end) = [];          % How many clicks did each ear actually hear?
Preliminary_Data.average_clicks(:,point+1:end) = [];      % Average number of clicks to be the mean of the poisson distribution
Preliminary_Data.number_of_frames = 60 * Preliminary_Data.stimulus_duration; % How many 'frames' or chances for a click per second in a trial?

Preliminary_Data.order_of_clicks(point+1:end,:,:) = [];

%% Save final data to folder
if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
    mkdir(fullfile(directory, 'RawData'));
    fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
    save(fileName, 'Preliminary_Data'); % save the data
else
    fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
    save(fileName, 'Preliminary_Data'); % save the data
end

end