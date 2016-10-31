function [] = ExperimentAuditory(subjectID, automatic, phase, add_noise, directory)

% Example Input - Experiment('Matthew', 'Matthew', 0, 0, '/Users/bcs206/Documents/Summer/')

if ~exist('automatic','var')  || ~exist('add_noise','var') || ~exist('directory','var')
    automatic = 0;       % 0 = normal subject, 1 = ideal observer auto running the experiment, 2 = human_like like observer, 3 = human observer
    phase = 0;           % 0 = volume experiment, 1 = ratio experiment
    add_noise = 0;       % 0 = no noise added, 1 = some background noise added to the sound
    directory = pwd;     % Make it equal to the current directory
end

% This functions initializes and runs an experiment using PsychToolBox

% subjectID is a string to dictate which subject is currently running the
% experiment. Ex. Experiment('01')

% automatic determines if it's to be a person or the computer running the experiment
% Ex. Experiment('01', 0) = person, Experiment('01', 1) = computer auto-runs as the ideal observer

% phase determines if the subject is being tested on volume/contrast or
% ratio. 0 = volume/contrast, 1 = ratio


% add_noise includes some background noise
% Ex. Experiment('01', 0, 0) = no noise added, Experiment('01', 0, 1) = some noise added

% directory allows this code to be able to create and save files of the subject data on any computer

settings = LoadSettings(directory);


%% Set Up the Initialization of the expeirment
cd(fullfile(directory, 'Code')) % Set the current directory
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

% Set up keyboard functions
KbName('UnifyKeyNames');
goKey = KbName(settings.keyGo);
exitKey = KbName(settings.keyExit);

% This is the first preliminary phase with a constant ratio (20, 4) and finding the threshold volume

if phase == 0
    
    fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        %% Instruction Screen
        Screen('TextSize', wPtr, 20); % Set text size to 20
        Screen('DrawText', wPtr, 'You will hear several clicks in each trial.', xc-500, yc-150, white);
        Screen('DrawText', wPtr, 'Then you will be asked which ear heard the higher rate of clicks.', xc-500, yc-100, white);
        Screen('DrawText', wPtr, 'Select the arrow key corresponding to the answer within 1 sec after trial ends.', xc-500, yc-50, white);
        Screen('DrawText', wPtr, 'A high pitched beep means correct, a low pitched beep means incorrect.', xc-500, yc, white);
        Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+50, white);
        Screen('DrawText', wPtr, sprintf('Press %s to begin.', settings.keyGoName), xc-500, yc+100, white);
        Screen('Flip', wPtr); % Function to flip to the next screen image
        [~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
        while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
            [~, ~, keyCode] = KbCheck;
        end
        Screen('Flip', wPtr); % Function to flip to the next screen image
        
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials =100;
        loops = 2;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.volume = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);  % Re ords underlying click rate
        Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);%stores answers based on actual click rates for staircase
        
        Preliminary_Data.current_trial = 0;
        Preliminary_Data.test_phase = ([1:loops].*preliminary_trials) - preliminary_trials + 1;
        
        % Set properties of a click
        Preliminary_Data.sampling_rate = 6000;   % Indirectly controls frequency
        Preliminary_Data.stimulus_duration = 2;           % How many seconds is each trial?
        Preliminary_Data.bins = Preliminary_Data.sampling_rate * Preliminary_Data.stimulus_duration;         % There needs to be an number of bins for the sound vector
        
        Preliminary_Data.click_rate = zeros(2,preliminary_trials*loops);          % How many clicks did each ear actually hear?
        Preliminary_Data.average_clicks = zeros(2,preliminary_trials*loops);      % Average number of clicks to be the mean of the poisson distribution
        Preliminary_Data.number_of_frames = 60 * Preliminary_Data.stimulus_duration; % How many 'frames' or chances for a click per second in a trial?
        
        Preliminary_Data.order_of_clicks = zeros(preliminary_trials*loops, 2, Preliminary_Data.number_of_frames);
        % Store all clicks sounded, ever
        % Number of trials and the sampling rate times stimulus duration to store the two sounds for the two ears
        
        max_volume = 0.5;      % Starting background volume level
        volume = max_volume;
        move_on = 0;        % Did the subject get two correct trials yet?
        step_size = 1.5;    % How strongly should the volume level be adjusted?
        previous_trial = 1;     % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks how many reversals the subject has gotten so far
        
        % Begin Preliminary Trials
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
                        Screen('DrawText', wPtr, sprintf('Press %s whenever you are ready again.', settings.keyGoName), xc-500, yc-50, white);
                        
                        Screen('Flip', wPtr); % Function to flip to the next screen image
                        [~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
                        while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
                            [~, ~, keyCode] = KbCheck;
                        end
                        Screen('Flip', wPtr);
                    end
                end
                volume = max_volume;
                move_on = 0;
                step_size = 1.5;
                previous_trial = 1;
                reversal_counter = 0;
                %max_volume = 0.1;
            else
                flag=0;
            end
            
            Preliminary_Data.current_trial = i;
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.volume(i) = volume;
            
            if rand < 0.5   % The left ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = 20;
                Preliminary_Data.average_clicks(2,i) = 4;
            else            % The right ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = 4;
                Preliminary_Data.average_clicks(2,i) = 20;
            end
            
            Preliminary_Data.ratio(i) = 20;
            
            % The correct anxwer is based on underlying click rate, not on the actual number of clicks each ear hears
            if Preliminary_Data.average_clicks(1,i) > Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif Preliminary_Data.average_clicks(1,i) < Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.correct_answer(i) = 1;
                else
                    Preliminary_Data.correct_answer(i) = 0;
                end
            end
            
            % Pass in the screen being used, subject ID, the struct with
            % all of the data, the trial numbers, and the fact it's the person or computer
            % running the experiment
            I = trialStimuliAuditory(wPtr, subjectID, Preliminary_Data, i, automatic, phase, add_noise, directory, settings);
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice;     % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.click_rate(1,i) = I.number_of_left_clicks;
            Preliminary_Data.click_rate(2,i) = I.number_of_right_clicks;   % Actual # of clicks for each frame
            
            Preliminary_Data.order_of_clicks(i, :, :) = I.clicks; % Record binary vector indicating the prescence or lack of click in each time bin
            
            % The staircase is based on the actual click rate, not on the underlying number of clicks each ear hears
            if Preliminary_Data.click_rate(1,i) > Preliminary_Data.click_rate(2,i)
                Preliminary_Data.staircase_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif Preliminary_Data.click_rate(1,i) < Preliminary_Data.click_rate(2,i)
                Preliminary_Data.staircase_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.staircase_answer(i) = 1;
                else
                    Preliminary_Data.staircase_answer(i) = 0;
                end
            end
            %% Feedback & Accuracy
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 0)
                Preliminary_Data.accuracy(i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1, 0.2);              % Beep for correct when it's the person running the experiment
                end
                
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 1)
                Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
            elseif (isnan(Preliminary_Data.choice(i)) )
                %Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(2, 0.2);              % Buzz for 'invalid trial'
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
                pause(.6); % Pause for 600 ms after feedback before next trial
            end
            if reversal_counter == 8       % It's a rule of thumb for when it's time to change the step size, since the subject has reversed 5 times
                reversal_counter = 0;  % Reset counter
                if step_size == 1.5
                    step_size = 1.2;
                elseif step_size == 1.2
                    step_size = 1.1;
                else
                    Preliminary_Data.test_phase(ceil(i/preliminary_trials)) = i; % This is when the preliminary phase ends and the test phase data starts
                    reversal_counter = 9;
                end
            end
            
            if move_on == 0 && ~isnan(Preliminary_Data.choice(i))
                if volume < max_volume
                    volume = volume*step_size;		% Subject got the trial wrong and the volume needs to be increased
                    %volume = volume/step_size; % as step sixe less than 1
                else
                    volume = max_volume;  % Just an insurance measure to keep the volume from going too high
                end
            elseif move_on == 2 && ~isnan(Preliminary_Data.choice(i))
                move_on = 0;  
                %volume = volume*step_size;  % As step size less than 1
                volume = volume/step_size;   % Subject got two trials right and it's time to decrease the volume level
            end
            % if move_on is equal to 1, then nothing needs to be changed
            
            %% Save the data after every trial in case of shut-downs
            if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
                mkdir(fullfile(directory, 'RawData'));
                fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            else
                fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            end
            if ~isnan(Preliminary_Data.choice(i))
                i=i+1;
            end
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        end
        %{
        test_start = 1;
        test_end = 0;
        Test_Data.current_trial = 0;
        for t = 1:loops
            stage = Preliminary_Data.test_phase(t);
            test_end =  test_end + preliminary_trials*t - stage + 1;
            
            Test_Data.move_on(test_start:test_end) = Preliminary_Data.move_on(stage:preliminary_trials*t);
            Test_Data.step_size(test_start:test_end) = Preliminary_Data.step_size(stage:preliminary_trials*t);
            Test_Data.reversal_counter(test_start:test_end) = Preliminary_Data.reversal_counter(stage:preliminary_trials*t);
            Test_Data.volume(test_start:test_end) = Preliminary_Data.volume(stage:preliminary_trials*t);
            Test_Data.correct_answer(test_start:test_end) = Preliminary_Data.correct_answer(stage:preliminary_trials*t);
            Test_Data.reaction_time(test_start:test_end) = Preliminary_Data.reaction_time(stage:preliminary_trials*t);
            Test_Data.choice(test_start:test_end) = Preliminary_Data.choice(stage:preliminary_trials*t);
            Test_Data.accuracy(test_start:test_end) = Preliminary_Data.accuracy(stage:preliminary_trials*t);
            Test_Data.ratio(test_start:test_end) = Preliminary_Data.ratio(stage:preliminary_trials*t);
            Test_Data.current_trial = Test_Data.current_trial + (preliminary_trials*t - stage + 1);
            Test_Data.test_phase = Preliminary_Data.test_phase;
            Test_Data.sampling_rate = Preliminary_Data.sampling_rate;
            Test_Data.stimulus_duration = Preliminary_Data.stimulus_duration;
            Test_Data.bins = Preliminary_Data.bins;
            Test_Data.click_rate(:,test_start:test_end) = Preliminary_Data.click_rate(:,stage:preliminary_trials*t);
            Test_Data.average_clicks(:,test_start:test_end) = Preliminary_Data.average_clicks(:,stage:preliminary_trials*t);
            Test_Data.number_of_frames = Preliminary_Data.number_of_frames;
            Test_Data.order_of_clicks(test_start:test_end,:,:) = Preliminary_Data.order_of_clicks(stage:preliminary_trials*t,:,:);
            Test_Data.staircase_answer(test_start:test_end) = Preliminary_Data.staircase_answer(stage:preliminary_trials*t);
            
            test_start = test_end+1;
        
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryTestVolume.mat']); % create a name for the data you want to save
            save(fileName, 'Test_Data'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryTestVolume.mat']); % create a name for the data you want to save
            save(fileName, 'Test_Data'); % save the data
        end
        %}
    end
    
elseif phase == 1
    
    % This is the second preliminary phase with a constant volume (max volume of 1) and finding the threshold ratio
    
    fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        %% Instruction Screen
        Screen('TextSize', wPtr, 20); % Set text size to 20
        Screen('DrawText', wPtr, 'You will hear several clicks in each trial.', xc-500, yc-150, white);
        Screen('DrawText', wPtr, 'Then you will be asked which ear heard the higher rate of clicks.', xc-500, yc-100, white);
        Screen('DrawText', wPtr, 'Select the arrow key corresponding to the answer within 2 secs after trial ends.', xc-500, yc-50, white);
        Screen('DrawText', wPtr, 'A high pitched beep means correct, a low pitched beep means incorrect.', xc-500, yc, white);
        Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+50, white);
        Screen('DrawText', wPtr, sprintf('Press %s to begin.', settings.keyGoName), xc-500, yc+100, white);
        Screen('Flip', wPtr); % Function to flip to the next screen image
        [~, ~, keyCode] = KbCheck;      % Variable to track the next keyboard press
        while ~keyCode(goKey)        % While loop to wait for the spacebar to be pressed
            [~, ~, keyCode] = KbCheck;
        end
        Screen('Flip', wPtr); % Function to flip to the next screen image
        
        %% Preliminary Calibration Phase
        
        % Set up struct to store data/answers
        preliminary_trials = 100;
        loops = 4;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.volume = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);  % Re ords underlying click rate
        Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);
        
        Preliminary_Data.current_trial = 0;
        Preliminary_Data.test_phase = ([1:loops].*preliminary_trials) - preliminary_trials + 1;
        
        % Set properties of a click
        Preliminary_Data.sampling_rate = 6000;   % Indirectly controls frequency
        Preliminary_Data.stimulus_duration = 2;           % How many seconds is each trial?
        Preliminary_Data.bins = Preliminary_Data.sampling_rate * Preliminary_Data.stimulus_duration;         % There needs to be an number of bins for the sound vector
        
        Preliminary_Data.click_rate = zeros(2,preliminary_trials*loops);          % How many clicks did each ear actually hear?
        Preliminary_Data.average_clicks = zeros(2,preliminary_trials*loops);      % Average number of clicks to be the mean of the poisson distribution
        Preliminary_Data.number_of_frames = 60 * Preliminary_Data.stimulus_duration; % How many 'frames' or chances for a click per second in a trial?
        
        Preliminary_Data.order_of_clicks = zeros(preliminary_trials*loops, 2, Preliminary_Data.number_of_frames);
        % Store all clicks sounded, ever
        % Number of trials and the sampling rate times stimulus duration to store the two sounds for the two ears
        
        ratio_sum = 24;
        max_ratio = 20;      % Starting ratio of clicks, 20 for one ear and 4 for the other ear
        ratio = max_ratio;
        step_size = 1;    % How strongly should the volume level be adjusted?
        move_on = 0;        % Did the subject get two correct trials yet?
        previous_trial = 1;     % Tracks if the subject was right or wrong last trial (1 is right and 0 is wrong on previous_trial trial)
        reversal_counter = 0;	% Tracks how many reversals the subject has gotten so far
        
        % Begin Preliminary Trials
        flag=0;
        i=1;
        while i <= preliminary_trials * loops
            
            if mod(i,preliminary_trials) == 1
                if i~=1
                    flag=flag+1;
                    if automatic == 0 && flag==1
                        sounds(-1, 1.5);
                        Screen('TextSize', wPtr, 20); % Set text size to 20
                        Screen('DrawText', wPtr, 'You finished a block.', xc-500, yc-150, white);
                        Screen('DrawText', wPtr, 'You may take a break!', xc-500, yc-100, white);
                        Screen('DrawText', wPtr, sprintf('Press %s whenever you are ready again.', settings.keyGoName), xc-500, yc-50, white);
                        
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
                step_size = 1;
                previous_trial = 1;
                reversal_counter = 0;
            else
                flag=0;
            end
            
            Preliminary_Data.current_trial = i;
            
            Preliminary_Data.move_on(i) = move_on;
            Preliminary_Data.step_size(i) = step_size;
            Preliminary_Data.volume(i) = 1;
            Preliminary_Data.reversal_counter(i) = reversal_counter;
            Preliminary_Data.ratio(i) = ratio;
            
            if rand < 0.5   % The left ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = ratio;
                Preliminary_Data.average_clicks(2,i) = ratio_sum - ratio;
            else            % The right ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = ratio_sum - ratio;
                Preliminary_Data.average_clicks(2,i) = ratio;
            end
            % The correct anxwer is based on underlying click rate, not on the actual number of clicks each ear hears
            if Preliminary_Data.average_clicks(1,i) > Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif Preliminary_Data.average_clicks(1,i) < Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.correct_answer(i) = 1;
                else
                    Preliminary_Data.correct_answer(i) = 0;
                end
            end
            
            % Pass in the screen being used, subject ID, the struct with
            % all of the data, the trial numbers, and the fact it's the person or computer
            % running the experiment
            I = trialStimuliAuditory(wPtr, subjectID, Preliminary_Data, i, automatic, phase, add_noise, directory, settings);
            
            % The correct anxwer is based on underlying click rate, not on the actual number of clicks each ear hears
            if Preliminary_Data.average_clicks(1,i) > Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif Preliminary_Data.average_clicks(1,i) < Preliminary_Data.average_clicks(2,i)
                Preliminary_Data.correct_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.correct_answer(i) = 1;
                else
                    Preliminary_Data.correct_answer(i) = 0;
                end
            end
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice;     % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.click_rate(1,i) = I.number_of_left_clicks;
            Preliminary_Data.click_rate(2,i) = I.number_of_right_clicks;   % Actual # of clicks for each frame
            
            Preliminary_Data.order_of_clicks(i, :, :) = I.clicks; % Record binary vector indicating the prescence or lack of click in each time bin
            % The staircase is based on the actual click rate, not on the underlying number of clicks each ear hears
            if Preliminary_Data.click_rate(1,i) > Preliminary_Data.click_rate(2,i)
                Preliminary_Data.staircase_answer(i) = 1;     % If 1, the answer was left, and if 0, the answer was right
            elseif Preliminary_Data.click_rate(1,i) < Preliminary_Data.click_rate(2,i)
                Preliminary_Data.staircase_answer(i) = 0;
            else % the ears have an equal underlying click rate
                if rand < 0.5
                    Preliminary_Data.staircase_answer(i) = 1;
                else
                    Preliminary_Data.staircase_answer(i) = 0;
                end
            end
            %% Feedback & Accuracy
            if (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 1) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 0)
                Preliminary_Data.accuracy(i) = 1;	% 1 is true for accuracy
                if automatic == 0
                    sounds(1, 0.2);              % Beep for correct when it's the person running the experiment
                end
                
            elseif (Preliminary_Data.choice(i) == 1 && Preliminary_Data.correct_answer(i) == 0) || (Preliminary_Data.choice(i) == 0 && Preliminary_Data.correct_answer(i) == 1)
                Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(0, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
            elseif (isnan(Preliminary_Data.choice(i)) )
                %Preliminary_Data.accuracy(i) = 0;	% 0 is false for inaccuracy
                if automatic == 0
                    sounds(2, 0.2);              % Buzz for wrong when it's the person running the experiment
                end
            end
            %% Staircase & Accuracy
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
            
            %% Save the data after every trial in case of shut-downs
            if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
                mkdir(fullfile(directory, 'RawData'));
                fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            else
                fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            end
            if ~isnan(Preliminary_Data.choice(i))
                i=i+1;
            end
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        end
        %{
        test_start = 1;
        test_end = 0;
        Test_Data.current_trial = 0;
        for t = 1:loops
            stage = Preliminary_Data.test_phase(t);
            test_end =  test_end + preliminary_trials*t - stage + 1;
            
            Test_Data.move_on(test_start:test_end) = Preliminary_Data.move_on(stage:preliminary_trials*t);
            Test_Data.step_size(test_start:test_end) = Preliminary_Data.step_size(stage:preliminary_trials*t);
            Test_Data.reversal_counter(test_start:test_end) = Preliminary_Data.reversal_counter(stage:preliminary_trials*t);
            Test_Data.volume(test_start:test_end) = Preliminary_Data.volume(stage:preliminary_trials*t);
            Test_Data.correct_answer(test_start:test_end) = Preliminary_Data.correct_answer(stage:preliminary_trials*t);
            Test_Data.reaction_time(test_start:test_end) = Preliminary_Data.reaction_time(stage:preliminary_trials*t);
            Test_Data.choice(test_start:test_end) = Preliminary_Data.choice(stage:preliminary_trials*t);
            Test_Data.accuracy(test_start:test_end) = Preliminary_Data.accuracy(stage:preliminary_trials*t);
            Test_Data.ratio(test_start:test_end) = Preliminary_Data.ratio(stage:preliminary_trials*t);
            Test_Data.current_trial = Test_Data.current_trial + (preliminary_trials*t - stage + 1);
            Test_Data.test_phase = Preliminary_Data.test_phase;
            Test_Data.sampling_rate = Preliminary_Data.sampling_rate;
            Test_Data.stimulus_duration = Preliminary_Data.stimulus_duration;
            Test_Data.bins = Preliminary_Data.bins;
            Test_Data.click_rate(:,test_start:test_end) = Preliminary_Data.click_rate(:,stage:preliminary_trials*t);
            Test_Data.average_clicks(:,test_start:test_end) = Preliminary_Data.average_clicks(:,stage:preliminary_trials*t);
            Test_Data.number_of_frames = Preliminary_Data.number_of_frames;
            Test_Data.order_of_clicks(test_start:test_end,:,:) = Preliminary_Data.order_of_clicks(stage:preliminary_trials*t,:,:);
            Test_Data.staircase_answer(test_start:test_end) = Preliminary_Data.staircase_answer(stage:preliminary_trials*t);
            
            test_start = test_end+1;
        
        end
        
        %% Save final data to folder
        if ~exist(fullfile(directory, 'RawData'), 'dir') % Check the directory actually exists
            mkdir(fullfile(directory, 'RawData'));
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryTestRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Test_Data'); % save the data
        else
            fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryTestRatio.mat']); % create a name for the data you want to save
            save(fileName, 'Test_Data'); % save the data
        end
        %}
    end
end

Screen('CloseAll'); % close screen
end