function [] = ExperimentGabor(subjectID_prelim, subjectID, automatic, directory)

% Example Input - ExperimentGabor('Matthew', 0, '/Users/bcs206/Documents/Summer/')

if ~exist('automatic','var') || ~exist('directory','var')
    automatic = 0;       % 0 = normal subject, 1 = ideal observer auto running the experiment
    directory = pwd;     % Make it equal to the current directory
end

% This functions initializes and runs an experiment using PsychToolBox

% subjectID is a string to dictate which subject is currently running the
% experiment. Ex. Experiment('01')

% automatic determines if it's to be a person or the computer running the experiment
% Ex. Experiment('01', 0) = person, Experiment('01', 1) = computer auto-runs the idela observer

% directory allows this code to be able to create and save files of the subject data on any computer


%% Set Up the Initialization of the expeirment
cd(fullfile(directory, 'Code')) % Set the current directory
directory = fullfile(directory, 'RawData');  % Directory to save the data and files to
commandwindow; % Moves the cursor to the commandwindow

InitializeMatlabOpenGL


% Screen set up
whichScreen = 0; %allow to choose the display if there's more than one
ResolutionScreen = Screen('Resolution', whichScreen); % Gets screen resolution
ScreenSize = [0 0 ResolutionScreen.width ResolutionScreen.height]; % Sets full screen
xc = ScreenSize(3)/2; %	Gets the middle of the horizontal axis
yc = ScreenSize(4)/2; % Gets the middle of the vertical axis
Screen('Preference', 'SkipSyncTests', 0); % Opens Screen

white = [255 255 255];          % Sets the color to be white
black = [0 0 0];                % Sets the color to be black

load('ColorCal2MeasurementsNoHotMirror.mat');
gammaTable = gammaTable2;
% Screen('LoadNormalizedGammaTable', win, gammaTable*[1 1 1]);
[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [], 32); % Opens window, sets background as black, sets screensize
%[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [xc-900 yc-500 xc+900 yc+500], 32);
% Creates a small window instead of using the full screen
% Mainly to allow screenshots

Screen('LoadNormalizedGammaTable', wPtr, gammaTable*[1 1 1]);

% Set up keyboard functions
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
left = KbName('leftArrow');
right = KbName('rightArrow');
up = KbName('upArrow');
down = KbName('downArrow');

if automatic == 0     % If automatic == 1, skip the preliminary phase since it's the computer running the experiment
    
    % This is the first preliminary phase with a constant ratio (10, 0) and finding the threshold contrast between 255 and 127
	
	fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryContrast.mat']); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
		Screen('DrawText', wPtr, 'You will see a series of images flashing very quickly in the middle of the screen.', xc-500, yc-150, white);
		Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target below the images.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
		Screen('DrawText', wPtr, sprintf('Select the image positioned to the left or right by pressing %s or %s respectively', settings.keyLeftName, settings.keyRightName), xc-500, yc+50, white);
		Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
		Screen('DrawText', wPtr, 'Press the space bar to begin.', xc-500, yc+150, white);    % Display text colored white
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
		
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 20;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.contrast = zeros(1,preliminary_trials);         % How clear an image is, or the signal level
        Preliminary_Data.number_of_images = 10;                                 % How many images/frames shown in a trial
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials);         % What was the right orientation/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials);                 % Which orientation did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials);                % Did the subject make the right orientation selection? 1 = yes, 0 = no
        Preliminary_Data.order_of_orientations = zeros(preliminary_trials,Preliminary_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
        Preliminary_Data.log_odds = zeros(1,preliminary_trials);
        Preliminary_Data.ratio = zeros(1,preliminary_trials);
        
        max_contrast = 128.0;     % Starting contrast level
        contrast = max_contrast;
        move_on = 0;        % Did the subject get two correct trials yet?
        step_size = 2.0;    % How strongly should the contrast level be adjusted?
        previous_trial = 1;           % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks of how many reversals the subject has gotten so far
        % It does this by counting when previous_trial switches from 0 to 1 or 1 to 0
        
        %{
        % Set properties of a gabor patch
        Preliminary_Data.gabor_sigma_x = 0.6;	% how large the patch is along the x-axis
        Preliminary_Data.gabor_sigma_y = 0.6;	% how large the patch is along the y-axis
        Preliminary_Data.gabor_spatial_frequency = 5;            % the thickness of the bands
        Preliminary_Data.gabor_phase = 0;         % phase
        Preliminary_Data.gabor_startpoint = -1;   % image size
        Preliminary_Data.gabor_step = 20;
        Preliminary_Data.gabor_endpoint = 1;
        %}
        % Properties which vary from trial to trial
        Preliminary_Data.screen_frame = 6;	% how long each image will show on screen in frame rates
        Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
        Preliminary_Data.image_length_x = 4;  % Size of the image along x-axis
        Preliminary_Data.image_length_y = 4;
        
        Preliminary_Data.image_template1 = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template2 = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template_difference = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        
        %image_collection = zeros(preliminary_trials,Preliminary_Data.number_of_images,Preliminary_Data.gabor_step,Preliminary_Data.gabor_step);
        image_collection = zeros(preliminary_trials, Preliminary_Data.number_of_images, ...
            Preliminary_Data.image_length_x*Preliminary_Data.screen_resolution, Preliminary_Data.image_length_y*Preliminary_Data.screen_resolution);
        
        
        % Begin the first phase of Preliminary Trials
        for i = 1:preliminary_trials
            
            Preliminary_Data.current_trial = i;
            % In other functions, it'll be useful to know the trial number without excess inputs
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.contrast(i) = contrast;
            
            
            % This sets the orientations, the probabilities for each
            % orientation, and their ratios to each orientation
            if rand <= 0.5
                desired_orientations = 1;        % Left Orientation
                probabilities = [1];             % 100% chance of only one orientation
                Preliminary_Data.ratio(i) = 1;  % Use only the one orientation
            else
                desired_orientations = 0;        % Right Orientation
                probabilities = [1];             % 100% chance of only one orientation
                Preliminary_Data.ratio(i) = 1;  % Use only the one orientation
            end
            
            % Function to generate the orientations of all images in the
            % trial and determine correct orientation to select
            [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Preliminary_Data.number_of_images);
            
            % A function to generate a struct containing properties of the images used in this trial, I,
            % and an collection of the gabors which will be displayed to the screen, image_array
            image_array = makeImages(Preliminary_Data, order_of_orientations, contrast);
            
            % Store all images shown
            image_collection(i,:,:,:) = image_array;
            
            % Pass in the contrast level, all of the images, screen being
            % used, subject ID, the struct with all of the data, and the
            % fact it's the person or computer running the experiment
            I = trialStimuliGabor(contrast, image_array, wPtr, subjectID, Preliminary_Data, automatic, directory);
            
            Preliminary_Data.correct_answer(i) = correct_orientation_answer;
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.log_odds(i) = I.log_odds;
            
            for k = 1:Preliminary_Data.number_of_images
                Preliminary_Data.image_template1(i,k) = I.log_regress(1,k);   % Log odds for left orientation image
                Preliminary_Data.image_template2(i,k) = I.log_regress(2,k);   % Log odds for right orientation image
                Preliminary_Data.image_template_difference(i,k) = I.log_regress(3,k);   % Log odds for difference in above two images image
            end
            
            
            %% Feedback & Accuracy
            if (I.choice == 1 && I.log_odds > 0) || (I.choice == 0 && I.log_odds < 0)
                Preliminary_Data.accuracy(1,i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1);              % Beep for correct when it's the person running the experiment
                end
                move_on = move_on + 1;	% if right, we increment move_on and later check it to decrease contrast level
                if previous_trial == 0    % if the subject got the last trial wrong
                    previous_trial = 1;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            elseif (I.choice == 0 && I.log_odds > 0) || (I.choice == 1 && I.log_odds < 0)
                Preliminary_Data.accuracy(1,i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0);              % Buzz for wrong when it's the person running the experiment
                end
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the contrast level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            end
            
            Preliminary_Data.order_of_orientations(i,:) = order_of_orientations;  % Record random ordering of all orientations
            
            if reversal_counter == 15       % It's a rule of thumb for when it's time to change the step size, since the subject has reversed 15 times
                reversal_counter = 0; % Reset counter
                if step_size == 2.0
                    step_size = 1.5;
                elseif step_size == 1.5
                    step_size = 1.2;
                elseif step_size == 1.2      % Note that we never go back up a step size
                    step_size = 1.1;
                else
                    break;           % Subjects have reached threshold
                end
            end
            
            if move_on == 0
                if contrast < max_contrast
                    contrast = contrast*step_size;		% Subject got the trial wrong and the contrast needs to be increased
                else
                    contrast = max_contrast;
                end
            elseif move_on == 2
                move_on = 0;                    	% Subject got two trials right and it's time to lower the contrast level
                contrast = contrast/step_size;
            end
            % if move_on is equal to 1, then nothing needs to be changed
            
            if contrast > max_contrast
                contrast = max_contrast;  % Just an insurance measure to keep the contrast from going over the maximum value allowed
            end
			
			% Maybe save the data after every nth trial? Remember, you're saving images too, so there's a noticiable delay between trials for the subjects
			
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryContrast.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryContrast.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        end
    end
	
    % This is the second preliminary phase with a constant contrast (223) and finding the threshold ratio between (10, 0) and (5, 5)
	
	fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryRatio.mat']); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
        Screen('DrawText', wPtr, 'You are now halfway through the preliminary phase.', xc-500, yc-200, white);
        Screen('DrawText', wPtr, 'This time the image brightness won''t change, but the image will switch between orientations.', xc-500, yc-150, white);
        Screen('DrawText', wPtr, 'You will see several images in each trial.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
		Screen('DrawText', wPtr, sprintf('Select the image positioned to the left or right by pressing %s or %s respectively', settings.keyLeftName, settings.keyRightName), xc-500, yc+50, white);
		Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
		Screen('DrawText', wPtr, 'Press the space bar to begin.', xc-500, yc+150, white);    % Display text colored white
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
		
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 25;
		block_trials = 20;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*block_trials);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*block_trials);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.contrast = zeros(1,preliminary_trials*block_trials);         % How clear an image is, or the signal level
        Preliminary_Data.number_of_images = 10;                                 % How many images/frames shown in a trial
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*block_trials);         % What was the right orientation/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*block_trials);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*block_trials);                 % Which orientation did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*block_trials);                % Did the subject make the right orientation selection? 1 = yes, 0 = no
        Preliminary_Data.order_of_orientations = zeros(preliminary_trials*block_trials,Preliminary_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
        Preliminary_Data.log_odds = zeros(1,preliminary_trials*block_trials);
        Preliminary_Data.ratio = zeros(1,preliminary_trials*block_trials);
        
		Preliminary_Data.current_trial = 0;
        Preliminary_Data.current_block = 0;
		
        ratio_sum = 10;
        max_ratio = 10;      % Starting ratio of orientations, 10 for one image and 0 for the other
        ratio = max_ratio;
        step_size = 1;    		% How strongly should the ratio be adjusted?
        
        % Properties which vary from trial to trial
        Preliminary_Data.screen_frame = 6;	% how long each image will show on screen in frame rates
        Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
        Preliminary_Data.image_length_x = 4;  % Size of the image along x-axis
        Preliminary_Data.image_length_y = 4;
        
        Preliminary_Data.image_template1 = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template2 = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        Preliminary_Data.image_template_difference = zeros(preliminary_trials,Preliminary_Data.number_of_images);
        
        %image_collection = zeros(preliminary_trials,Preliminary_Data.number_of_images,Preliminary_Data.gabor_step,Preliminary_Data.gabor_step);
        image_collection = zeros(preliminary_trials, Preliminary_Data.number_of_images, ...
            Preliminary_Data.image_length_x*Preliminary_Data.screen_resolution, Preliminary_Data.image_length_y*Preliminary_Data.screen_resolution);
        
        
        % Begin Preliminary Trials
        for i = 1:preliminary_trials
			
			Preliminary_Data.current_block = i;
            
            for block = 1:block_trials
				
                Preliminary_Data.current_trial = (i*20-20) + block;
                Preliminary_Data.ratio((i*20-20) + block) = ratio/max_ratio;		% Only records the % for the image that appears more frequently
                
                Preliminary_Data.step_size((i*20-20) + block) = step_size;
                Preliminary_Data.contrast((i*20-20) + block) = 223;   % The constant contrast level throughout this preliminary phase
                
				% This sets the orientations, the probabilities for each
				% orientation, and their ratios to each orientation
				if rand <= 0.5
					desired_orientations = [1, 0];        % Left Orientation, Right Orientation (Ordering is important!)
					probabilities = [ratio/max_ratio, (max_ratio - ratio)/max_ratio];		% Calculates % of images for one orientation vs another orientation
					Preliminary_Data.correct_answer = 1;	% Since left appears >= than the right orientation, it's the correct answer
				else
					desired_orientations = [0, 1];        % Right Orientation, Left Orientation (Ordering is important!)
					probabilities = [ratio/max_ratio, (max_ratio - ratio)/max_ratio];		% Calculates % of images for one orientation vs another orientation
					Preliminary_Data.correct_answer = 0;	% Since right appears >= than the left orientation, it's the correct answer
				end
				
				% Function to generate the orientations of all images in the
				% trial and determine correct orientation to select
				[order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Preliminary_Data.number_of_images);
				
				% A function to generate a struct containing properties of the images used in this trial, I,
				% and an collection of the gabors which will be displayed to the screen, image_array
				image_array = makeImages(Preliminary_Data, order_of_orientations, contrast);
				
				% Store all images shown
				image_collection(i,:,:,:) = image_array;
				
				% Pass in the contrast level, all of the images, screen being
				% used, subject ID, the struct with all of the data, and the
				% fact it's the person or computer running the experiment
				I = trialStimuliGabor(contrast, image_array, wPtr, subjectID, Preliminary_Data, automatic, directory);
				
				Preliminary_Data.correct_answer((i*20-20) + block) = correct_orientation_answer;
				Preliminary_Data.reaction_time((i*20-20) + block) = I.reaction;
				Preliminary_Data.choice((i*20-20) + block) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
				Preliminary_Data.log_odds((i*20-20) + block) = I.log_odds;
				
				for k = 1:Preliminary_Data.number_of_images
					Preliminary_Data.image_template1((i*20-20) + block,k) = I.log_regress(1,k);   % Log odds for left orientation image
					Preliminary_Data.image_template2((i*20-20) + block,k) = I.log_regress(2,k);   % Log odds for right orientation image
					Preliminary_Data.image_template_difference((i*20-20) + block,k) = I.log_regress(3,k);   % Log odds for difference in above two images image
				end
				
				%% Feedback & Accuracy
				if (Preliminary_Data.choice((i*20-20) + block) == 1 && Preliminary_Data.correct_answer((i*20-20) + block) == 1) || (Preliminary_Data.choice((i*20-20) + block) == 0 && Preliminary_Data.correct_answer((i*20-20) + block) == 0)
					Preliminary_Data.accuracy((i*20-20) + block) = 1;	% 1 is true for accuracy
					if automatic == 0
						sounds(1);              % Beep for correct when it's the person running the experiment
					end
				elseif (Preliminary_Data.choice((i*20-20) + block) == 0 && Preliminary_Data.correct_answer((i*20-20) + block) == 1) || (Preliminary_Data.choice((i*20-20) + block) == 1 && Preliminary_Data.correct_answer((i*20-20) + block) == 0)
					Preliminary_Data.accuracy((i*20-20) + block) = 0;	% 0 is false for inaccuracy
					if automatic == 0
						sounds(0);              % Buzz for wrong when it's the person running the experiment
					end
				end
				
				Preliminary_Data.order_of_orientations(i,:) = order_of_orientations;  % Record random ordering of all orientations
            end
            
            accuracy_rating = sum(Preliminary_Data.accuracy((i*20-19):(i*20)))/20;
            
            if accuracy_rating >= 0.8           % Subject did well and it's time to lower the ratio
                ratio = ratio - step_size;
                if ratio < ratio_sum/2
                    ratio = ratio_sum/2;  % To set a minimum level for the staircase
                end
            elseif accuracy_rating <= 0.7       % Subject did poorly and it's time to increase the ratio
                if ratio < max_ratio
                    ratio = ratio + step_size;
                end
                if ratio > max_ratio
                    ratio = max_ratio;  % Just an insurance measure to keep the ratio from going too high
                end
            end
            % if accuracy_rating is 0.7-0.8, then nothing needs to be changed
			
			% Maybe save the data after every nth trial? Remember, you're saving images too, so there's a noticiable delay between trials for the subjects
			
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
        end
    end
    % Otherwise skip the entire prelim right to loading the prelims file for the needed threshold volume and ratios
    
    % Get threshold contrast
    fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryContrast.mat']);
    results = load(fileName);
    % This will save to a struct with the field, results.Preliminary_Data
    
    contrast_graph = results.Preliminary_Data.contrast;
    
    average_last_trials = 10;
    average_contrast = 0;
    
    for i = 0:average_last_trials      % Add last 10 trials of the prelims and divide by 10
        average_contrast = average_contrast + contrast_graph(i + (results.Preliminary_Data.current_trial - average_last_trials));
    end
    contrast = average_contrast/average_last_trials; % threshold contrast in test phase
    
    % Get threshold ratio
    fileName = fullfile(directory, 'RawData', [subjectID '-GaborPreliminaryRatio.mat']);
    results = load(fileName);
    % This will save to a struct with the field, results.Preliminary_Data
    
    ratio_graph = results.Preliminary_Data.ratio;
    
    average_last_trials = 10;
    average_ratio = 0;
    
    for i = 0:average_last_trials      % Add last 10 trials of the prelims and divide by 10
        average_ratio = average_ratio + ratio_graph(i + (results.Preliminary_Data.current_trial - average_last_trials));
    end
    ratio = average_volume/average_last_trials;
    ratio = round(ratio);                       % threshold ratio in test phase
else
    contrast = 1.0;   % When running hte experiment for the computer, just set the contrast to an abitrary value
	ratio = 6;
end

%% Instruction Screen
Screen('TextSize', wPtr, 20); % Set text size to 20
Screen('DrawText', wPtr, 'You have completed the preliminary calibrations! Congratulations!', xc-500, yc-200, white);
Screen('DrawText', wPtr, 'You will see a series of images flashing very quickly in the middle of the screen.', xc-500, yc-150, white);
Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target below the images.', xc-500, yc-100, white);
Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
Screen('DrawText', wPtr, sprintf('Select the image positioned to the left or right by pressing %s or %s respectively', settings.keyLeftName, settings.keyRightName), xc-500, yc+50, white);
Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
Screen('DrawText', wPtr, 'Press the space bar to begin.', xc-500, yc+150, white);    % Display text colored white
Screen('Flip', wPtr); % Function to flip to the next screen image
[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
    [~, ~, keyCode] = KbCheck;
end
Screen('Flip', wPtr); % Function to flip to the next screen image

%% Test Phase

test_trials = 20;

Test_Data.contrast = zeros(1,test_trials);   % How clear an image is (0, 6]
Test_Data.number_of_images = 10;              % How many images shown in a trial
Test_Data.correct_answer = zeros(1,test_trials);   % What was the right orientation?
Test_Data.reaction_time = zeros(1,test_trials);    % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Test_Data.choice = zeros(1,test_trials);           % Which orientation did the subject choose? 1 = left, 0 = right
Test_Data.accuracy = zeros(1,test_trials);         % Did the subject make the right orientation selection? 1 = yes, 0 = no
Test_Data.order_of_orientations = zeros(test_trials,Test_Data.number_of_images);   % Record the randomized order of image orientations throughout the trial
Test_Data.log_odds = zeros(1,test_trials);
Test_Data.ratio = zeros(1,test_trials);

Test_Data.image_length_x = 4;
Test_Data.image_length_y = 4;

% Properties which vary from trial to trial
Test_Data.screen_frame = 6;	% how long each image will show on screen in frame rates
Test_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor

Test_Data.image_template1 = zeros(test_trials,Test_Data.number_of_images);
Test_Data.image_template2 = zeros(test_trials,Test_Data.number_of_images);
Test_Data.image_template_difference = zeros(test_trials,Test_Data.number_of_images);

image_collection = zeros(test_trials,Test_Data.number_of_images,Test_Data.screen_resolution*4,Test_Data.screen_resolution*4);

fileName = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']);
if exist(fileName, 'file')
    subjectID = sprintf('%sI',subjectID);
    load(fileName);
    contrast = min(Test_Data.contrast);
end

% Begin Test Trials
for i = 1:test_trials

    if i == test_trials/2 +1
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
		Screen('DrawText', wPtr, 'You are now halfway through the experiment.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'While the task will be identical as before,', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'the brightness of the images will be much lower and harder to see.', xc-500, yc, white);
		Screen('DrawText', wPtr, 'Press the spacebar to continue.', xc-500, yc+50, white);
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
	end
    
    Test_Data.current_trial = i;
    
    if i <= test_trials/2		% High contrast (223) and threshold ratio condition
		Test_Data.contrast = 223;
		
		ratio_sum = Test_Data.number_of_images;
		if rand < 0.5
			desired_orientations = [1, 0];   % [Left, Right]
			probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];
		else
			desired_orientations = [0, 1];   % [Right, Left]
			probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];
		end
        Test_Data.ratio(i) = ratio/ratio_sum;
        [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Test_Data.number_of_images);
    else						% Threshold contrast and high ratio (10,0) condition
		Test_Data.contrast = contrast;
		if rand <= 0.5
            desired_orientations = [1];   % Left orientation image only
            probabilities = [1];
        else
            desired_orientations = [0];  % Right orientation image only
            probabilities = [1];
		end
        Test_Data.ratio(i) = ratio_sum;
        [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Test_Data.number_of_images);
    end
	
    image_array = makeImages(Test_Data, order_of_orientations, contrast);
    % A function to generate a struct containing properties of the images used in this trial, I,
    % and an collection of the gabors which will be displayed to the screen, image_array
    
    image_collection(i,:,:,:) = image_array;
    % Store all images shown, ever
    
    I = trialStimuliGabor(contrast, image_array, wPtr, subjectID, Test_Data, automatic, directory);
    
    Test_Data.correct_answer(i) = correct_orientation_answer;
    Test_Data.reaction_time(i) = I.reaction;
    Test_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
    Test_Data.log_odds(i) = I.log_odds;
    
    %% Feedback & Accuracy
    if ((I.choice == 1 && I.log_odds > 0)  ||  (I.choice == 0 && I.log_odds < 0))
        % Column 13: Accuracy of Orientation Choice
        Test_Data.accuracy(1,i) = 1;	% 1 is true for accuracy
        if automatic == 0   % only sound when there's a subject
            sounds(1);              % Beep for correct
        end
    else % if ((I.choice == 0 && I.log_odds > 0)  ||  (I.choice == 1 && I.log_odds < 0))
        % Column 13: Accuracy of Orientation Choice
        Test_Data.accuracy(1,i) = 0;	% 0 is false for inaccuracy
        if automatic == 0   % only sound when there's a subject
            sounds(0);              % Buzz for wrong
        end
    end
    
    % Column 14+: Orientations
    Test_Data.order_of_orientations(i,:) = order_of_orientations;
    
    % Column : Logistic Regressions
    for k = 1:Test_Data.number_of_images
        Test_Data.image_template1(i,k) = I.log_regress(1,k);
        Test_Data.image_template2(i,k) = I.log_regress(2,k);
        Test_Data.image_template_difference(i,k) = I.log_regress(3,k);
    end
    
	%% Save the data after every 25th trial in case of shut-downs
	if mod(i, 25) == 0
		if ~exist(directory, 'dir')
			mkdir(directory);

			fileName = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']); % create a name for the data you want to save as a csv
			save(fileName, 'Test_Data', 'image_collection'); % save the data
		else
			fileName = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']); % create a name for the data you want to save as a csv
			save(fileName, 'Test_Data', 'image_collection'); % save the data
		end
	end
end


%% Save final data to folder
if ~exist(fullfile(directory, 'RawData'), 'dir')
    mkdir(fullfile(directory, 'RawData'));
    fileName = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']); % create a name for the data you want to save as a csv
    save(fileName, 'Test_Data', 'image_collection'); % save the data
else
    fileName = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']); % create a name for the data you want to save as a csv
    save(fileName, 'Test_Data', 'image_collection'); % save the data
end

Screen('CloseAll'); % close screen
end