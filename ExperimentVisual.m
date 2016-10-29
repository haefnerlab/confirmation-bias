function [] = ExperimentVisual(subjectID, automatic, directory)

% Example Input - ExperimentVisual('Matthew', 0, '/Users/bcs206/Documents/Summer/')

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
cd([directory 'Code/']) % Set the current directory
directory = [directory 'RawData/'];  % Directory to save the data and files to
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

[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [], 32); % Opens window, sets background as black, sets screensize
%[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [xc-900 yc-500 xc+900 yc+500], 32);
% Creates a small window instead of using the full screen
% Mainly to allow screenshots

% Set up keyboard functions
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
%% Instruction Screen
left = KbName('leftArrow');
right = KbName('rightArrow');
up = KbName('upArrow');
down = KbName('downArrow');


Screen('TextSize', wPtr, 20); % Set text size to 20
Screen('DrawText', wPtr, 'You will see two bars flashing very quickly in the middle of the screen with a background of static.', xc-500, yc-150, white);
Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target between the bars', xc-500, yc-100, white);
Screen('DrawText', wPtr, 'Then you will be asked which side had the more frequently appearing bar.', xc-500, yc-50, white);
Screen('DrawText', wPtr, 'Select the left or right side by pressing the corresponding arrow key.', xc-500, yc, white);
Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+50, white);
Screen('DrawText', wPtr, 'Press the spacebar to begin.', xc-500, yc+100, white);    % Display text colored white
Screen('Flip', wPtr); % Function to flip to the next screen image
[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
    [~, ~, keyCode] = KbCheck;
end
Screen('Flip', wPtr); % Function to flip to the next screen image

if automatic == 0     % If automatic == 1, skip the preliminary phase since it's the computer running the experiment
    fileName = sprintf('%s%s-VisualPreliminary.mat',[directory 'RawData/'],subjectID); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 20;
        
        max_contrast = 127.0;    % Starting contrast level (background gray/128 + 127 = white/255)
        contrast = max_contrast;
        move_on = 0;        % Did the subject get two correct trials yet?
        step_size = 2.0;    % How strongly should the contrast level be adjusted?
        previous_trial = 1;     % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks how many reversals the subject has gotten so far
        % It does this by counting when previous_tri
        Preliminary_Data.number_of_images = 1;                       % How many images/frames shown in a trial
        Preliminary_Data.move_on = zeros(1,preliminary_trials);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.contrast = zeros(1,preliminary_trials);         % How bright a bar is, or the signal level
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials);         % What was the right bar/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials);                 % Which bar did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.order_of_flashes = zeros(preliminary_trials,2,Preliminary_Data.number_of_images);
                                                   % Record the randomized order of flashes throughout the trial in the premlims, it only appears once
        Preliminary_Data.screen_frame = 1;        % How long will each image show on the screen
        Preliminary_Data.image_length_x = 1;      % The width of the bar in ratio to the height
        Preliminary_Data.image_length_y = 10;     % The height of the bar in ration to the width
        Preliminary_Data.screen_resolution = 50;   % Magnify the size of the bar, or how many pixels to one datapoint of the image
        
        
        % Begin Preliminary Trials
        for i = 1:preliminary_trials
            
            Preliminary_Data.current_trial = i;
            % In other functions, it'll be useful to know the trial number without excess inputs
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.contrast(i) = contrast;
            
            % Set the Gaussian noise ahead of time
            left_noise = randn*16;
            right_noise = randn*16;
            
            % Decide if the bar will flash left or right in the bar detection task
            % There will only be one image fo the prelims and 20 for the test phase
            if rand < 0.5
                L = ones(1,Preliminary_Data.number_of_images).*(127 + contrast + left_noise);   % The bar will flash a color between gray and white (background + contrast + noise level)
                R = ones(1,Preliminary_Data.number_of_images).*(127 + right_noise);                % The bar will flash background gray + noise
                Preliminary_Data.correct_answer(i) = 1;   % Record correct answer
            else
                L = ones(1,Preliminary_Data.number_of_images).*(127 + left_noise);
                R = ones(1,Preliminary_Data.number_of_images).*(127 + contrast + right_noise);
                Preliminary_Data.correct_answer(i) = 0;
            end
            
            % If any of the images were greater than the maximum value allowed, set them to 255
            for im = 1:Preliminary_Data.number_of_images
                if L(im) > 255, L(im) = 255; end
                if R(im) > 255, R(im) = 255; end
            end
            
            % Record the sequence of bar color values
            Preliminary_Data.order_of_flashes(i,1,:) = L';
            Preliminary_Data.order_of_flashes(i,2,:) = R';
            
            % Pass in the screen being used, subject ID, the struct with
            % all of the data, and the fact it's the person or computer
            % running the experiment
            I = trialStimuliVisual(wPtr, subjectID, Preliminary_Data, automatic, directory);
            
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            
            Preliminary_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
            
            
           %% Feedback & Accuracy
            if ((I.choice == 1 && Preliminary_Data.correct_answer(i) == 1)  ||  (I.choice == 0 && Preliminary_Data.correct_answer(i) == 0))
                Preliminary_Data.accuracy(i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1, 0.2);              % Beep for correct when it's the person running the experiment
                end
                move_on = move_on + 1;	% If right, we increment move_on and later check it to decrease contrast level
                if previous_trial == 0    % If the subject got the last trial wrong
                    previous_trial = 1;       % then the staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            elseif ((I.choice == 1 && Preliminary_Data.correct_answer(i) == 0)  ||  (I.choice == 0 && Preliminary_Data.correct_answer(i) == 1))
                Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
                move_on = 0;            % if wrong, we reset move_on to 0 and later increase the contrast level
                if previous_trial == 1    % if the subject got the last trial right
                    previous_trial = 0;       % staircase is reversing
                    reversal_counter = reversal_counter+1;
                end
            end
            
            
            if reversal_counter == 15       % It's a rule of thumb for when it's time to change the step size, since the subject has reversed 15 times
                reversal_counter = 0;   % Reset counter
                if step_size == 2
                    step_size = 1.5;
                elseif step_size == 1.5
                    step_size = 1.2;
                elseif step_size == 1.2
                    step_size = 1.1;
                else
                    break;           % Subjects have reached threshold
                end
            end
            
            if move_on == 0
                if contrast < max_contrast
                    contrast = round(contrast*step_size);		% Subject got the trial wrong and the contrast needs to be lowered
                end
            elseif move_on == 2
                move_on = 0;                    	% Subject got two trials right and it's time to increase the contrast level
                contrast = round(contrast/step_size);
            end
            % if move_on is equal to 1, then nothing needs to be changed
            
            if contrast > max_contrast
                contrast = max_contrast;  % Just an insurance measure to keep the contrast from going over the maximum value allowed
            end
            
            %% Save the data after every trial
            if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
                mkdir([directory 'RawData/']);
                fileName = sprintf('%s%s-VisualPreliminary.mat',directory,subjectID); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            else
                fileName = sprintf('%s%s-VisualPreliminary.mat',directory,subjectID); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            end
            
        end
    end
    % Otherwise if the subject has already done the prelims, skip the entire prelim phase
    % right to loading the prelims file for the needed threshold contrast
    % level (average of last 10 trials)
    
    load(fileName);
    
    contrast_graph = Preliminary_Data.contrast;
    
    average_last_trials = 10;
    
    average_contrast = 0;
    for i = 1:average_last_trials      % Add last 10 trials of the prelims and divide by 10
        average_contrast = average_contrast + contrast_graph(1,i+(Preliminary_Data.current_trial-average_last_trials));
    end
    
    contrast = average_contrast/average_last_trials;
	
	fileName = sprintf('%s%s-VisualTest.mat',[directory 'RawData/'],subjectID);
	while exist(fileName, 'file')
		subjectID = [subjectID, 'I'];
		fileName = sprintf('%s%s-VisualTest.mat',[directory 'RawData/'],subjectID);
	end
else
    contrast = 1.0;   % When running the experiment for the computer, just set the contrast to an abitrary value
end

Screen('TextSize', wPtr, 20); % Set text size to 20
Screen('DrawText', wPtr, 'You have completed the preliminary calibrations.', xc-500, yc-250, white);
Screen('DrawText', wPtr, 'This time the brightness of the bars will not change over the course of multiple trials.', xc-500, yc-200, white);
Screen('DrawText', wPtr, 'You will see two bars flashing very quickly in the middle of the screen with a background of static.', xc-500, yc-150, white);
Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target between the bars', xc-500, yc-100, white);
Screen('DrawText', wPtr, 'Then you will be asked which side had the more frequently appearing bar.', xc-500, yc-50, white);
Screen('DrawText', wPtr, 'Select the left or right side by pressing the corresponding arrow key.', xc-500, yc, white);
Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+50, white);
Screen('DrawText', wPtr, 'Press the spacebar to begin.', xc-500, yc+100, white);    % Display text colored white
Screen('Flip', wPtr); % Function to flip to the next screen image
[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
    [~, ~, keyCode] = KbCheck;
end
Screen('Flip', wPtr); % Function to flip to the next screen image


%% Test Phase

test_trials = 40;

Test_Data.number_of_images = 30;    % How many images/frames shown in a trial (EAch iamge will be a two frame-pair)
Test_Data.average_flashes = zeros(2,test_trials);  % How many flahes should appear on average?
Test_Data.contrast = zeros(1,test_trials);         % How bright a bar is, or the signal level
Test_Data.correct_answer = zeros(1,test_trials);         % What was the right bar/answer?
Test_Data.reaction_time = zeros(1,test_trials);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Test_Data.choice = zeros(1,test_trials);                 % Which bar did the subject choose? 1 = left, 0 = right
Test_Data.accuracy = zeros(1,test_trials);               % Did the subject make the right selection? 1 = yes, 0 = no
Test_Data.order_of_flashes = zeros(test_trials,2, Test_Data.number_of_images);
                   % Record the randomized order of flashes throughout the trial in the for each bar
Test_Data.screen_frame = 1;        % How long will each image show on the screen
Test_Data.image_length_x = 1;      % The width of the bar in ratio to the height
Test_Data.image_length_y = 10;     % The height of the bar in ration to the width
Test_Data.screen_resolution = 50;   % Magnify the size of the bar, or how many pixels to one datapoint of the image
Test_Data.flash_rate = zeros(2,test_trials);   % How many times did the bars actually flash for each side


% Begin Test Trials
for i = 1:test_trials

    if i == (test_trials/2 +1)
		%% Instruction Screen
		Screen('TextSize', wPtr, 20); % Set text size to 20
		Screen('DrawText', wPtr, 'You are now halfway through the experiment.', xc-500, yc-100, white);
		Screen('DrawText', wPtr, 'While the task will be identical as before,', xc-500, yc-50, white);
		Screen('DrawText', wPtr, 'the brightness of the bar will be much lower and harder to see.', xc-500, yc, white);
		Screen('DrawText', wPtr, 'Press the spacebar to continue.', xc-500, yc+50, white);
		Screen('Flip', wPtr); % Function to flip to the next screen image
		[~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
		while ~keyCode(spaceKey)        % While loop to wait fo rhte spacebar to be pressed
			[~, ~, keyCode] = KbCheck;
		end
		Screen('Flip', wPtr); % Function to flip to the next screen image
	end
    
    Test_Data.current_trial = i;
    Test_Data.contrast(i) = contrast;
    
    L = ones(1,Test_Data.number_of_images*2);
    R = ones(1,Test_Data.number_of_images*2);
    
    if i <= test_trials/2     % For the first half of the trials, both bars will flash 10 times on average
        Test_Data.average_flashes(1,i) = 10;
        Test_Data.average_flashes(2,i) = 10;
    else                      % For the second half of the trials, one side will flash 20 times on average
        if rand < 0.5
            Test_Data.average_flashes(1,i) = 20;
            Test_Data.average_flashes(2,i) = 0;
        else
            Test_Data.average_flashes(1,i) = 0;
            Test_Data.average_flashes(2,i) = 20;
        end
    end
    
    flashes_L = 0; % Count nubmer of flashes for the left and right bar
    flashes_R = 0;
    
    % For each image, they will be given a color value
    for p = 1:Test_Data.number_of_images*2
        if i <= test_trials/2  % The high-contrast condition
            
            left_noise = 0;%randn*16;   % The mean is 0 with 16 sigma
            right_noise = 0;%randn*16;  % The same noise needs to be added to the frame-pair
            if mod(p,2) == 0       % If on the even trials, generate stimulus
                left_odds = Test_Data.average_flashes(1,i)/Test_Data.number_of_images;   % 10 flashes/30 screen frames = 1/3
                right_odds = Test_Data.average_flashes(2,i)/Test_Data.number_of_images;
                L(p) = rand < left_odds;       % The dice roll to determine if there will be a flash or not, 1 = yes 0 = no
                R(p) = rand < right_odds;
                if L(p) == 1    % Convert 1s and 0s to flashes
                    L(p) = (223 + left_noise);
                    flashes_L = flashes_L + 1;    % Increment flash counter
                else
                    L(p) = (127 + left_noise);
                end
                if R(p) == 1    % Convert 1s and 0s to flashes
                    R(p) = (223 + right_noise);
                    flashes_R = flashes_R + 1;    % Increment flash counter
                else
                    R(p) = (127 + right_noise);
                end
                if L(p) > 255, L(p) = 255; end  % If greater than max value, set to max value
                if R(p) > 255, R(p) = 255; end
            else
                L(p) = 127 + left_noise;   % Every other stimulus image will appear as background gray plus noise, or on the odd numbered trials
                R(p) = 127 + right_noise;
            end
            
        else                   % The low-contrast condition
            
            left_noise = randn*16;   % The mean is 0 with 16 sigma
            right_noise = randn*16;  % The same noise needs to be added to both frames in the frame-pair
            if mod(p,2) == 0
                left_odds = Test_Data.average_flashes(1,i)/Test_Data.number_of_images;   % 20 flashes/30 screen frames = 2/3
                right_odds = Test_Data.average_flashes(2,i)/Test_Data.number_of_images;  % 0 flashes/30 screen frames = 0/3
                L(p) = rand < left_odds;       % The dice roll to determine if there will be a flash or not, 1 = yes 0 = no
                R(p) = rand < right_odds;
                if L(p) == 1    % Convert 1s and 0s to flashes
                    L(p) = (127 + contrast + left_noise);
                    flashes_L = flashes_L + 1;    % Increment flash counter
                else
                    L(p) = (127 + left_noise);
                end
                if R(p) == 1
                    R(p) = (127 + contrast + right_noise);
                    flashes_R = flashes_R + 1;    % Increment flash counter
                else
                    R(p) = (127 + right_noise);
                end
                if L(p) > 255, L(p) = 255; end  % If greater than max value, set to max value
                if R(p) > 255, R(p) = 255; end
            else
                L(p) = 127 + left_noise;   % Every other stimulus image will be the background color plus noise
                R(p) = 127 + right_noise;
            end
            
        end
    end
    Test_Data.order_of_flashes(i,1,:) = L(2:2:Test_Data.number_of_images*2)';  % Don't bother save the stimuli matching the background
    Test_Data.order_of_flashes(i,2,:) = R(2:2:Test_Data.number_of_images*2)';  % Save only the even-numbered trials
    
    Test_Data.flash_rate(1,i) = flashes_L;  % Save the number of flashes for each side
    Test_Data.flash_rate(2,i) = flashes_R;
    
    if flashes_L == flashes_R   % If there are an equal number of flashes, randomly pick a sied
        if rand < 0.5
            Test_Data.correct_answer(i) = 1;
        else
            Test_Data.correct_answer(i) = 0;
        end
    else
        Test_Data.correct_answer(i) = (flashes_L > flashes_R);  % If 1, the left image flashes more. If 0, it's the right image
    end
    
    
    % Pass in the screen being used, subject ID, the struct with
    % all of the data, and the fact it's the person or computer
    % running the experiment
    I = trialStimuliVisual(wPtr, subjectID, Test_Data, automatic, directory);
    
    Test_Data.reaction_time(i) = I.reaction;
    
    Test_Data.choice(i) = I.choice; % If 1, subject chose left, and if 0, the subject chose right
    
    
    %% Feedback & Accuracy
    if ((I.choice == 1 && Test_Data.correct_answer(i) == 1)  ||  (I.choice == 0 && Test_Data.correct_answer(i) == 0))
        Test_Data.accuracy(i) = 1;	% 1 is true for accuracy
        if automatic == 0   % only sound when there's a subject
            sounds(1, 0.2);              % Beep for correct
        end
    elseif ((I.choice == 1 && Test_Data.correct_answer(i) == 0)  ||  (I.choice == 0 && Test_Data.correct_answer(i) == 1))
        Test_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
        if automatic == 0   % only sound when there's a subject
            sounds(0, 0.2);              % Buzz for wrong
        end
    end
    
    
    %% Save the data after every trial
    if ~exist([directory 'RawData/'], 'dir')
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-VisualTest.mat',directory,subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-VisualTest.mat',directory,subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
end
%%%% May not need this part hence commented!!!!!!!!!!!!!!!!!
%% Save final data to folder
%if ~exist([directory 'RawData/'], 'dir')
   % mkdir([directory 'RawData/']);
    
    %fileName = sprintf('%s%s-VisualTest.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
    %save(fileName, 'Test_Data'); % save the data
%else
    %fileName = sprintf('%s%s-VisualTest.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
    %save(fileName, 'Test_Data'); % save the data
%end

Screen('CloseAll'); % close screen
end