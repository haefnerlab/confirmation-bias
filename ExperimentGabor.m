function [] = ExperimentGabor(subjectID, automatic, phase, directory, varargin)

% Example Input - ExperimentGabor('Matthew', 0, '/Users/bcs206/Documents/Summer/')

if ~exist('automatic','var') || ~exist('directory','var')
    automatic = 0;       % 0 = normal subject, 1 = ideal observer auto running the experiment
    phase = 0;           % 0 = contrast experiment, 1 = ratio experiment
    directory = pwd;     % Make it equal to the current directory
end

% This functions initializes and runs an experiment using PsychToolBox

% subjectID is a string to dictate which subject is currently running the
% experiment. Ex. Experiment('01')

% automatic determines if it's to be a person or the computer running the experiment
% Ex. Experiment('01', 0) = person, Experiment('01', 1) = computer auto-runs the idela observer

% directory allows this code to be able to create and save files of the subject data on any computer

settings = LoadSettings(directory);

%% Set Up the Initialization of the expeirment
cd([directory 'Code\']) % Set the current directory
directory = [directory 'RawData\'];  % Directory to save the data and files to
commandwindow; % Moves the cursor to the commandwindow

if settings.useOpenGL, InitializeMatlabOpenGL; end

% Screen set up
whichScreen = settings.whichScreen; %allow to choose the display if there's more than one
ResolutionScreen = Screen('Resolution', whichScreen); % Gets screen resolution
ScreenSize = [0 0 ResolutionScreen.width ResolutionScreen.height]; % Sets full screen
xc = ScreenSize(3)/2; %	Gets the middle of the horizontal axis
yc = ScreenSize(4)/2; % Gets the middle of the vertical axis
Screen('Preference', 'SkipSyncTests', settings.ptbSkipSyncTests); % Opens Screen

white = [255 255 255];          % Sets the color to be white
black = [0 0 0];                % Sets the color to be black

[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [], 32); % Opens window, sets background as black, sets screensize
%[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [xc-900 yc-500 xc+900 yc+500], 32);
% Creates a small window instead of using the full screen
% Mainly to allow screenshots
Screen('LoadNormalizedGammaTable', wPtr, gammaTable*[1 1 1]);

if ~isempty(settings.gammaTableFile)
    gtdata = load(settings.gammaTableFile);
    Screen('LoadNormalizedGammaTable', wPtr, gtdata.(settings.gammaTable)*[1 1 1]);
end

tracker_info = EyeTracker.initEyeTracker(whichScreen, ...
    'fixationSymbol', 'b', ...
    'fixationCenter', [xc, ScreenSize(4)-50], ...
    'fixationSymbolSize', [30 30], ...
    varargin{:});

% Set up keyboard functions
KbName('UnifyKeyNames');
goKey = KbName(settings.keyGo);
exitKey = KbName(settings.keyExit);

% This is the first preliminary phase with a constant ratio (20, 4) and finding the threshold contrast
HideCursor(whichScreen)
if phase == 0
    
%     EyeTracker.AutoCalibrate(tracker_info);
    
    fileName = sprintf('%s%s-GaborDataContrast.mat',directory,subjectID); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
		Screen('DrawText', wPtr, 'You will see a series of images flashing very quickly in the middle of the screen.', xc-500, yc-150, white);
		Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target below the images.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
		Screen('DrawText', wPtr, 'Select the image positioned to the left or right by pressing the corresponding arrow key.', xc-500, yc+50, white);
		Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
		Screen('DrawText', wPtr, sprintf('Press %s to begin.', settings.goKeyName), xc-500, yc+150, white);    % Display text colored white
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
		
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 40;
        loops = 2;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.contrast = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
        Preliminary_Data.number_of_images = 10;
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.order_of_orientations = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
        Preliminary_Data.log_odds = zeros(1,preliminary_trials*loops);
        Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);
        Preliminary_Data.average_orientations = zeros(2,preliminary_trials*loops);
        
        Preliminary_Data.current_trial = 0;
        Preliminary_Data.test_phase = ([1:loops].*preliminary_trials) - preliminary_trials + 1;
        
        Preliminary_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
        Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
        Preliminary_Data.image_length_x = 5;  % Size of the image along x-axis
        Preliminary_Data.image_length_y = 5;
        
        Preliminary_Data.image_template1 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template2 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template_difference = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        
        Preliminary_Data.eye_tracker_points = {};
        
        res = Preliminary_Data.screen_resolution;
        Preliminary_Data.left_template = zeros(res * 5);
        for i=1:5
            Preliminary_Data.left_template((i-1)*res+1:i*res, (i-1)*res+1:i*res) = 1;
        end
        Preliminary_Data.right_template = rot90(Preliminary_Data.left_template);
        
        
        image_collection = zeros(preliminary_trials*loops, Preliminary_Data.number_of_images, ...
            Preliminary_Data.image_length_x*Preliminary_Data.screen_resolution, Preliminary_Data.image_length_y*Preliminary_Data.screen_resolution);
        
        max_contrast = 128.0;      % Starting background contrast level
        contrast = max_contrast;
        move_on = 0;        % Did the subject get two correct trials yet?
        step_size = 2.0;    % How strongly should the contrast level be adjusted?
        previous_trial = 1;     % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks how many reversals the subject has gotten so far
        
        % Begin Preliminary Trials
        EyeTracker.AutoCalibrate(tracker_info);
        flag=0;
        i=1;
        while i <= preliminary_trials * loops
            
            if mod(i,preliminary_trials) == 1
                if i~=1
                    flag=flag+1;
                    if automatic == 0 && flag==1
                        sounds(-1, 1.5);
                        Screen('TextSize', wPtr, 20); % Set text size to 20
                        Screen('DrawText', wPtr, 'You have completed a block.', xc-500, yc-150, white);
                        Screen('DrawText', wPtr, 'You may take a break if you want!', xc-500, yc-100, white);
                        Screen('DrawText', wPtr, sprintf('Press %s whenever you are ready again.', settings.goKeyName), xc-500, yc-50, white);
                        
                        Screen('Flip', wPtr); % Function to flip to the next screen image
                        [~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
                        while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
                            [~, ~, keyCode] = KbCheck;
                        end
                        Screen('Flip', wPtr);
                    end
                end
                contrast = max_contrast;
                move_on = 0;
                step_size = 2.0;
                previous_trial = 1;
                reversal_counter = 0;
            else
                flag=0;
            end
            
            Preliminary_Data.current_trial = i;
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.contrast(i) = contrast;
            
            % This sets the orientations, the probabilities for each
            % orientation, and their ratios to each orientation
            
            if rand < 0.5   % The left ear will hear more clicks
                Preliminary_Data.average_orientations(1,i) = Preliminary_Data.number_of_images;
                Preliminary_Data.average_orientations(2,i) = 0;
                desired_orientations = 1;        % Left Orientation
                probabilities = [1];             % 100% chance of only one orientation
                Preliminary_Data.correct_answer(i) = 1;	% Since left appears >= than the right orientation, it's the correct answer
                Preliminary_Data.staircase_answer(i) = 1;
                Preliminary_Data.ratio(i) = Preliminary_Data.number_of_images;  % Use only the one orientation
            else            % The right ear will hear more clicks
                Preliminary_Data.average_orientations(1,i) = 0;
                Preliminary_Data.average_orientations(2,i) = Preliminary_Data.number_of_images;
                desired_orientations = 0;        % Right Orientation
                probabilities = [1];             % 100% chance of only one orientation
                Preliminary_Data.correct_answer(i) = 0;	% Since right appears >= than the left orientation, it's the correct answer
                Preliminary_Data.staircase_answer(i) = 0;
                Preliminary_Data.ratio(i) = Preliminary_Data.number_of_images;  % Use only the one orientation
            end
            
            % Function to generate the orientations of all images in the
            % trial and determine correct orientation to select
            [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Preliminary_Data.number_of_images);
            
            % A function to generate a struct containing properties of the images used in this trial, I,
            % and an collection of the gabors which will be displayed to the screen, image_array
            image_array = makeImages(Preliminary_Data, order_of_orientations, contrast);
            
            % Store all images shown
            image_collection(i,:,:,:) = image_array;
            
            Preliminary_Data.order_of_orientations(i,:) = order_of_orientations;  % Record random ordering of all orientations
            
            % Pass in the contrast level, all of the images, screen being
            % used, subject ID, the struct with all of the data, and the
            % fact it's the person or computer running the experiment
            [I, eye_tracker_points, broke_fixation] = trialStimuliGabor(i, image_array, wPtr, subjectID, Preliminary_Data, automatic, phase, directory, tracker_info, settings);
            
            if broke_fixation
                Screen('Flip', wPtr);
                sounds(0, 0.2);
                pause(1);
                continue;
            end
            
            Preliminary_Data.eye_tracker_points{i} = eye_tracker_points;
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.log_odds(i) = I.log_odds;
            
            for k = 1:Preliminary_Data.number_of_images
                Preliminary_Data.image_template1(i,k) = I.log_regress(1,k);   % Log odds for left orientation image
                Preliminary_Data.image_template2(i,k) = I.log_regress(2,k);   % Log odds for right orientation image
                Preliminary_Data.image_template_difference(i,k) = I.log_regress(3,k);   % Log odds for difference in above two images image
            end
            
            % The staircase is based on the actual click rate, not on the underlying number of clicks each ear hears
            Preliminary_Data.staircase_answer(i) = correct_orientation_answer;     % If 1, the answer was left, and if 0, the answer was right
            
            %% Feedback & Accuracy
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 0)
                Preliminary_Data.accuracy(i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1, 0.2);              % Beep for correct when it's the person running the experiment
                end
                %{
                move_on = move_on + 1;	% if right, we increment move_on and later check it to decrease contrast level
                if previous_trial == 0    % if the subject got the last trial wrong
                    previous_trial = 1;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
                %}
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 1)
                Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
                %{
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the contrast level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
                %}
            elseif (isnan(Preliminary_Data.choice(i)) )
                %Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(2, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
            end
            
            %% Staircase method
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.staircase_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.staircase_answer(i) == 0)
                
                move_on = move_on + 1;	% if right, we increment move_on and later check it to decrease volume level
                if previous_trial == 0    % if the subject got the last trial wrong
                    previous_trial = 1;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.staircase_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.staircase_answer(i) == 1)
                
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the volume level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            end
            if automatic == 0
                pause(.5); % Pause for 500 ms after feedback before next trial
            end
            
            if reversal_counter == 15       % It's a rule of thumb for when it's time to change the step size, since the subject has reversed 15 times
                reversal_counter = 0; % Reset counter
                if step_size == 2.0
                    step_size = 1.5;
                elseif step_size == 1.5
                    step_size = 1.2;
                elseif step_size == 1.2      % Note that we never go back up a step size
                    step_size = 1.1;
                else
                    Preliminary_Data.test_phase(ceil(i/preliminary_trials)) = i; % This is when the preliminary phase ends and the test phase data starts
                    reversal_counter = 16;
                end
            end
            
            if move_on == 0 && ~isnan(Preliminary_Data.choice(i))
                if contrast < max_contrast
                    contrast = contrast*step_size;		% Subject got the trial wrong and the contrast needs to be increased
                else
                    contrast = max_contrast;
                end
            elseif move_on == 2 && ~isnan(Preliminary_Data.choice(i))
                move_on = 0;                    	% Subject got two trials right and it's time to lower the contrast level
                contrast = contrast/step_size;
            end
            % if move_on is equal to 1, then nothing needs to be changed
            
            if contrast > max_contrast
                contrast = max_contrast;  % Just an insurance measure to keep the contrast from going over the maximum value allowed
            end
            
            % Maybe save the data after every nth trial? Remember, you're saving images too, so there's a noticiable delay between trials for the subjects
            if ~isnan(Preliminary_Data.choice(i))
                i=i+1;
            end
        end
        
        %% Save final data to folder
        if ~exist(directory, 'dir') % Check the directory actually exists
            mkdir(directory);
            fileName = sprintf('%s%s-GaborDataContrast.mat',directory,subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        else
            fileName = sprintf('%s%s-GaborDataContrast.mat',directory,subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        end
    end
    
elseif phase == 1
    
    % This is the second preliminary phase with a constant contrast (max contrast of 223) and finding the threshold ratio
    
    fileName = sprintf('%s%s-GaborDataRatio.mat',directory,subjectID); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
		Screen('DrawText', wPtr, 'You will see a series of images flashing very quickly in the middle of the screen.', xc-500, yc-150, white);
		Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target below the images.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
		Screen('DrawText', wPtr, 'Select the image positioned to the left or right by pressing the corresponding arrow key.', xc-500, yc+50, white);
		Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
		Screen('DrawText', wPtr, sprintf('Press %s to begin.', settings.goKeyName), xc-500, yc+150, white);    % Display text colored white
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
		
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 200;
        loops = 4;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.contrast = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
        Preliminary_Data.number_of_images = 10;
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.order_of_orientations = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
        Preliminary_Data.log_odds = zeros(1,preliminary_trials*loops);
        Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);
        Preliminary_Data.average_orientations = zeros(2,preliminary_trials*loops);
        
        Preliminary_Data.current_trial = 0;
        Preliminary_Data.test_phase = ([1:loops].*preliminary_trials) - preliminary_trials + 1;
        
        Preliminary_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
        Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
        Preliminary_Data.image_length_x = 5;  % Size of the image along x-axis
        Preliminary_Data.image_length_y = 5;
        
        Preliminary_Data.image_template1 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template2 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template_difference = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
        
        %image_collection = zeros(preliminary_trials,Preliminary_Data.number_of_images,Preliminary_Data.gabor_step,Preliminary_Data.gabor_step);
        image_collection = zeros(preliminary_trials*loops, Preliminary_Data.number_of_images, ...
            Preliminary_Data.image_length_x*Preliminary_Data.screen_resolution, Preliminary_Data.image_length_y*Preliminary_Data.screen_resolution);
        
        ratio_sum = Preliminary_Data.number_of_images;
        max_ratio = Preliminary_Data.number_of_images;      % Starting ratio of clicks, 20 for one ear and 4 for the other ear
        ratio = max_ratio;
        step_size = 1;    % How strongly should the contrast level be adjusted?
        move_on = 0;        % Did the subject get two correct trials yet?
        previous_trial = 1;     % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks how many reversals the subject has gotten so far
        
        contrast = 64;      % High Contrast
        
        % Begin Preliminary Trials
        EyeTracker.AutoCalibrate(tracker_info);
        i=1;
        flag=0;
        while i <= preliminary_trials * loops
            
            if mod(i,preliminary_trials) == 1
                if i~=1
                    flag=flag+1;
                    if automatic == 0 && flag==1
                        sounds(-1, 1.5);
                        Screen('TextSize', wPtr, 20); % Set text size to 20
                        Screen('DrawText', wPtr, 'You finished a block.', xc-500, yc-150, white);
                        Screen('DrawText', wPtr, 'You may take a break!', xc-500, yc-100, white);
                        Screen('DrawText', wPtr, sprintf('Press %s whenever you are ready again.', settings.goKeyName), xc-500, yc-50, white);
                        
                        Screen('Flip', wPtr); % Function to flip to the next screen image
                        [~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
                        while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
                            [~, ~, keyCode] = KbCheck;
                        end
                        Screen('Flip', wPtr);
                    end
                end
                
                ratio = max_ratio;
                move_on = 0;
                step_size = 1.0;
                previous_trial = 1;
                reversal_counter = 0;
            end
            
            Preliminary_Data.current_trial = i;
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.contrast(i) = contrast;
            
            
            % This sets the orientations, the probabilities for each
            % orientation, and their ratios to each orientation
            
            if rand < 0.5   % The left ear will hear more clicks
                Preliminary_Data.average_orientations(1,i) = ratio;
                Preliminary_Data.average_orientations(2,i) = ratio_sum - ratio;
                desired_orientations = [1, 0];        % Left Orientation
                probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];		% Calculates % of images for one orientation vs another orientation
                Preliminary_Data.correct_answer(i) = 1;	% Since left appears >= than the right orientation, it's the correct answer
                Preliminary_Data.ratio(i) = ratio;  % Use only the one orientation
            else            % The right ear will hear more clicks
                Preliminary_Data.average_orientations(1,i) = ratio_sum - ratio;
                Preliminary_Data.average_orientations(2,i) = ratio;
                desired_orientations = [0, 1];        % Right Orientation
                probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];		% Calculates % of images for one orientation vs another orientation
                Preliminary_Data.correct_answer(i) = 0;	% Since right appears >= than the left orientation, it's the correct answer
                Preliminary_Data.ratio(i) = ratio;  % Use only the one orientation
            end
            
            % Function to generate the orientations of all images in the
            % trial and determine correct orientation to select
            [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Preliminary_Data.number_of_images);
            Preliminary_Data.correct_answer(i) = correct_orientation_answer;
            % A function to generate a struct containing properties of the images used in this trial, I,
            % and an collection of the gabors which will be displayed to the screen, image_array
            image_array = makeImages(Preliminary_Data, order_of_orientations, contrast);
            
            % Store all images shown
            image_collection(i,:,:,:) = image_array;
            
            Preliminary_Data.order_of_orientations(i,:) = order_of_orientations;  % Record random ordering of all orientations
            
            % Pass in the contrast level, all of the images, screen being
            % used, subject ID, the struct with all of the data, and the
            % fact it's the person or computer running the experiment
            I = trialStimuliGabor(i, image_array, wPtr, subjectID, Preliminary_Data, automatic, phase, directory, tracker_info, settings);
            
            % The staircase is based on the actual click rate, not on the underlying number of clicks each ear hears
            if sum(Preliminary_Data.order_of_orientations(i,:)) > Preliminary_Data.number_of_images/2
                Preliminary_Data.staircase_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif sum(Preliminary_Data.order_of_orientations(i,:)) < Preliminary_Data.number_of_images/2
                Preliminary_Data.staircase_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.staircase_answer(i) = 1;
                else
                    Preliminary_Data.staircase_answer(i) = 0;
                end
            end
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.log_odds(i) = I.log_odds;
            
            for k = 1:Preliminary_Data.number_of_images
                Preliminary_Data.image_template1(i,k) = I.log_regress(1,k);   % Log odds for left orientation image
                Preliminary_Data.image_template2(i,k) = I.log_regress(2,k);   % Log odds for right orientation image
                Preliminary_Data.image_template_difference(i,k) = I.log_regress(3,k);   % Log odds for difference in above two images image
            end
            
            %% Feedback & Accuracy
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 0)
                Preliminary_Data.accuracy(i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1, 0.2);              % Beep for correct when it's the person running the experiment
                end
                %{
                move_on = move_on + 1;	% if right, we increment move_on and later check it to decrease contrast level
                if previous_trial == 0    % if the subject got the last trial wrong
                    previous_trial = 1;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
                %}
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 1)
                Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
                %{
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the contrast level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
                %}
            elseif (isnan(Preliminary_Data.choice(i)) )
                %Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(2, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
            end
            
            %% Staircase and Accuracy
            
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.staircase_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.staircase_answer(i) == 0)
                
                move_on = move_on + 1;	% if right, we increment move_on and later check it to decrease volume level
                if previous_trial == 0    % if the subject got the last trial wrong
                    previous_trial = 1;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.staircase_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.staircase_answer(i) == 1)
                
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the volume level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            end
            if automatic == 0
                pause(.6); % Pause for 600 ms after feedback before next trial
            end
            if reversal_counter == 5       % It's a rule of thumb for when it's time to change the step size, since the subject has reversed 5 times
                reversal_counter = 0;  % Reset counter
                if step_size == 1
                    step_size = 0.5;
                elseif step_size == 0.5
                    step_size = 0.25;
                else
                    Preliminary_Data.test_phase(ceil(i/preliminary_trials)) = i; % This is when the preliminary phase ends and the test phase data starts
                    reversal_counter = 6;
                end
            end
            
            if move_on == 0 && ~isnan(Preliminary_Data.choice(i))
                if ratio < max_ratio
                    ratio = ratio + step_size;
                end
                if ratio > max_ratio
                    ratio = max_ratio;  % Just an insurance measure to keep the ratio from going too high
                end
            elseif move_on == 2 && ~isnan(Preliminary_Data.choice(i))
                move_on = 0;                    	% Subject got two trials right and it's time to decrease the volume level
                ratio = ratio - step_size;
                if ratio < ratio_sum/2
                    ratio = ratio_sum/2;  % To set a minimum level for the staircase
                end
            end
            % if move_on is equal to 1, then nothing needs to be changed
            if ratio > max_ratio
                ratio = max_ratio;  % Just an insurance measure to keep the contrast from going over the maximum value allowed
            end
			
            % Maybe save the data after every nth trial? Remember, you're saving images too, so there's a noticiable delay between trials for the subjects
            if ~isnan(Preliminary_Data.choice(i))
                i=i+1;
            end
        end
        
        %% Save final data to folder
        if ~exist(directory, 'dir') % Check the directory actually exists
            mkdir(directory);
            fileName = sprintf('%s%s-GaborDataRatio.mat',directory,subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        else
            fileName = sprintf('%s%s-GaborDataRatio.mat',directory,subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        end
    end
end
ShowCursor([],whichScreen)
Screen('CloseAll'); % close screen
end