function [] = SaveQuitDataGabor(subjectID, phase, loops_consider, preliminary_trials, directory)

if phase == 0       % load the volume data
    filename = fullfile(directory, 'RawData', [subjectID '-GaborDataContrast.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
elseif phase == 1       % load the ratio data
    filename = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']);
    if ~exist(filename, 'file')
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
        return;
    else
        load(filename); % Load Preliminary_Data
    end
end

point = loops_consider*preliminary_trials;
Preliminary_Data.move_on (point+1:end) = [];          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
Preliminary_Data.step_size (point+1:end) = [];        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
Preliminary_Data.reversal_counter (point+1:end) = [];   % How many trials has the subject got wrong? When to change the step size?
Preliminary_Data.contrast(point+1:end) = [];         % How loud the sound is, or the signal level
Preliminary_Data.number_of_images = 10;
Preliminary_Data.correct_answer(point+1:end) = [];         % What was the right ear/answer?
Preliminary_Data.staircase_answer(point+1:end) = [];         % What was the right ear/answer?
Preliminary_Data.reaction_time(point+1:end) = [];          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Preliminary_Data.choice(point+1:end) = [];                 % Which ear did the subject choose? 1 = left, 0 = right
Preliminary_Data.accuracy(point+1:end) = [];               % Did the subject make the right selection? 1 = yes, 0 = no
Preliminary_Data.order_of_orientations((point+1:end),:) = [];        % Record the randomized order of image orientations throughout the trial
Preliminary_Data.log_odds(point+1:end) = [];
Preliminary_Data.ratio(point+1:end) = [];
Preliminary_Data.average_orientations(:,point+1) = [];

Preliminary_Data.current_trial = 0;
Preliminary_Data.test_phase = ([1:loops].*preliminary_trials) - preliminary_trials + 1;

Preliminary_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
Preliminary_Data.image_length_x = 4;  % Size of the image along x-axis
Preliminary_Data.image_length_y = 4;

Preliminary_Data.image_template1 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
Preliminary_Data.image_template2 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
Preliminary_Data.image_template_difference = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);


image_collection = zeros(preliminary_trials*loops, Preliminary_Data.number_of_images, ...
    Preliminary_Data.image_length_x*Preliminary_Data.screen_resolution, Preliminary_Data.image_length_y*Preliminary_Data.screen_resolution);
end