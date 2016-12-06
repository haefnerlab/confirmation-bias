function [GaborData, image_collection] = ExperimentGabor(subjectID, GaborData, varargin)

directory = fullfile(pwd, '..');
settings = LoadSettings(directory);

datadir = fullfile(directory, 'RawData');
if ~exist(datadir, 'dir'), mkdir(datadir); end

%% Environment and PsychToolBox Initialization
cd(fullfile(directory, 'Code')) % Set the current directory
commandwindow; % Moves the cursor to the commandwindow

% if settings.useOpenGL, InitializeMatlabOpenGL; end

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
HideCursor(whichScreen);

if ~isempty(settings.gammaTableFile)
    gtdata = load(settings.gammaTableFile);
    Screen('LoadNormalizedGammaTable', wPtr, gtdata.(settings.gammaTable)*[1 1 1]);
end

% Set up eye tracker, passing through any remaining varargin
tracker_info = EyeTracker.initEyeTracker(whichScreen, ...
    'fixationSymbol', 'b', ...
    'fixationCenter', [xc, 50], ...
    'fixationSymbolSize', [30 30], ...
    'fixationRadius', 18, ...
    varargin{:});

% Set up keyboard functions
KbName('UnifyKeyNames');
goKey = KbName(settings.keyGo);
exitKey = KbName(settings.keyExit);

if isequal(GaborData.stair_fn, @GaborStaircase.contrast)
    fileName = fullfile(datadir, [subjectID '-GaborDataContrast.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataContrastQuit.mat']);
elseif isequal(GaborData.stair_fn, @GaborStaircase.ratio)
    fileName = fullfile(datadir, [subjectID '-GaborDataRatio.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataRatioQuit.mat']);
elseif isequal(GaborData.stair_fn, @GaborStaircase.pixel_noise)
    fileName = fullfile(datadir, [subjectID '-GaborDataNoise.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataNoiseQuit.mat']);
else
    warning('No staircase-specific suffix on save name!');
    fileName = fullfile(datadir, [subjectID '.mat']);
    fileNameQuit = fullfile(datadir, [subjectID 'Quit.mat']);
end
if exist(fileName, 'file')
    Screen('CloseAll');
    error('Data for %s already exists', fileName);
end

%% Begin experiment
EyeTracker.AutoCalibrate(tracker_info);

% Instruction Screen
Screen('TextSize', wPtr, 20); % Set text size to 20
Screen('FillRect', wPtr, 127);
Screen('DrawText', wPtr, 'You will see a series of images flashing very quickly at the bottom of the screen.', xc-500, yc-150, white);
Screen('DrawText', wPtr, 'You are required to keep your eyes on the bull''s eye target at the top of the screen.', xc-500, yc-100, white);
Screen('DrawText', wPtr, 'Then you will be shown two images.', xc-500, yc-50, white);
Screen('DrawText', wPtr, 'You will have to decide which image appeared more frequently.', xc-500, yc, white);
Screen('DrawText', wPtr, sprintf('Select the image positioned to the left or right by pressing %s or %s respectively', settings.keyLeftName, settings.keyRightName), xc-500, yc+50, white);
Screen('DrawText', wPtr, 'Ask the researcher if you need further clarification.', xc-500, yc+100, white);
Screen('DrawText', wPtr, sprintf('Press %s to begin.', settings.keyGoName), xc-500, yc+150, white);    % Display text colored white
Screen('Flip', wPtr); % Function to flip to the next screen image
if ptbWaitKey([goKey exitKey]) == exitKey
    Screen('CloseAll');
    return;
end
Screen('Flip', wPtr); % Function to flip to the next screen image

% Preallocate images
image_collection = zeros(GaborData.blocks * GaborData.trials_per_block, ...
    GaborData.number_of_images, ...
    GaborData.image_length_x, ...
    GaborData.image_length_y);

    function earlyQuit
        save(fileNameQuit, 'GaborData', 'image_collection');
        ShowCursor();
        Screen('CloseAll');
    end

% Begin Preliminary Trials
seen_block_notification = false;
trial = 1;
block_trial = 1;
% Using 'while' rather than 'for' since invalid trials (broke fixation or
% didn't respond in time) don't increment 'trial'.
try
    while trial <= GaborData.trials_per_block * GaborData.blocks

        %% Bookkeeping to set up trial

        GaborData.current_trial = trial;
        % Reset params at the start of each block
        if mod(trial, GaborData.trials_per_block) == 1
            % Display a message if this is the beginning of the second or
            % higher block.
            if trial ~= 1 && isempty(GaborData.model_observer) && ~seen_block_notification
                seen_block_notification = true;
                sounds(-1, 1.5);
                Screen('TextSize', wPtr, 20); % Set text size to 20
                Screen('DrawText', wPtr, 'You have completed a block.', xc-500, yc-150, white);
                Screen('DrawText', wPtr, 'You may take a break if you want!', xc-500, yc-100, white);
                Screen('DrawText', wPtr, sprintf('Press %s whenever you are ready again.', settings.keyGoName), xc-500, yc-50, white);
                Screen('Flip', wPtr);

                if ptbWaitKey([goKey exitKey]) == exitKey
                    earlyQuit;
                    return;
                end

                Screen('Flip', wPtr);
            end

            % Start of a block - set params to initial values.
            GaborData.streak(trial) = 0;
            GaborData.reversal_counter(trial) = 0;
            GaborData.contrast(trial) = GaborData.contrast(1);
            GaborData.ratio(trial) = GaborData.ratio(1);
            GaborData.pixel_noise(trial) = GaborData.pixel_noise(1);
            GaborData.step_size(trial) = GaborData.step_size(1);
            block_trial = 1;
        else
            seen_block_notification = false;

            GaborData.contrast(trial) = GaborData.contrast(trial-1);
            GaborData.ratio(trial) = GaborData.ratio(trial-1);
            GaborData.pixel_noise(trial) = GaborData.pixel_noise(trial-1);
            GaborData.step_size(trial) = GaborData.step_size(trial-1);

            % Count correct streak (with respect to the ideal observer's
            % answer, not the underlying distribution)
            if GaborData.ideal_answer(trial-1) == GaborData.choice(trial-1)
                GaborData.streak(trial) = GaborData.streak(trial-1) + 1;
            else
                GaborData.streak(trial) = 0;
            end

            % Count reversals
            if block_trial > 2 && sign(GaborData.streak(trial-1)) ~= sign(GaborData.streak(trial))
                GaborData.reversal_counter(trial) = GaborData.reversal_counter(trial-1) + 1;
            else
                GaborData.reversal_counter(trial) = GaborData.reversal_counter(trial-1);
            end

            % Every 10 reversals, halve the step size
            if GaborData.reversal_counter(trial) > 1 && mod(GaborData.reversal_counter(trial), 10) == 0
                GaborData.step_size(trial) = GaborData.step_size(trial-1) / 2;
            end

            % Apply the staircase
            GaborData = GaborData.stair_fn(GaborData);
        end

        %% Run this trial

        % Randomly set each frame match (or mismatch) the correct choice for
        % this trail, using the current 'ratio' to decide.
        match_frames = 1 * (rand(1, GaborData.number_of_images) <= GaborData.ratio(trial));

        % Choose whether correct answer this trial will be Left or Right
        if rand < 0.5
            GaborData.correct_answer(trial) = 1; % this trial is Left
            GaborData.average_orientations(1, trial) = GaborData.number_of_images * GaborData.ratio(trial);
            GaborData.average_orientations(2, trial) = GaborData.number_of_images * (1 - GaborData.ratio(trial));
            GaborData.order_of_orientations(trial, :) = match_frames;
        else
            GaborData.correct_answer(trial) = 0; % this trial is Right
            GaborData.average_orientations(1, trial) = GaborData.number_of_images * (1 - GaborData.ratio(trial));
            GaborData.average_orientations(2, trial) = GaborData.number_of_images * GaborData.ratio(trial);
            GaborData.order_of_orientations(trial, :) = 1 - match_frames;
        end

        % A function to generate a struct containing properties of the images used in this trial, I,
        % and an collection of the gabors which will be displayed to the screen, image_array
        image_array = makeImages(GaborData);

        % Store all images shown
        image_collection(trial,:,:,:) = image_array;

        % Calculate log odds at each frame, both for the category of that frame
        % independent of the prior, and for the decision (including the prior)
        [log_frame_odds, log_decision_odds] = GaborLogOdds(image_array, ...
            GaborData.left_template, GaborData.right_template, ...
            GaborData.contrast(trial), GaborData.pixel_noise(trial)^2, ...
            GaborData.ratio(trial));
        GaborData.log_frame_odds(trial, :) = log_frame_odds;
        GaborData.log_decision_odds(trial, :) = log_decision_odds;

        % Record answer of the ideal observer
        GaborData.ideal_answer(trial) = 1 * (sum(GaborData.log_decision_odds(trial, :)) > 0);

        if isempty(GaborData.model_observer)
            % Pass in the contrast level, all of the images, screen being
            % used, subject ID, the struct with all of the data, and the
            % fact it's the person or computer running the experiment
            [I, eye_tracker_points, broke_fixation, quit] = trialStimuliGabor(GaborData, image_array, wPtr, tracker_info, settings);

            if quit, earlyQuit; return; end

            if broke_fixation || isnan(I.choice)
                Screen('FillRect', wPtr, 127);
                Screen('Flip', wPtr);
                sounds(2, 0.2);
                pause(1);
                continue;
            end

            GaborData.choice(trial) = I.choice;
            GaborData.reaction_time(trial) = I.reaction;
            GaborData.eye_tracker_points{trial} = eye_tracker_points;
        elseif strcmpi(GaborData.model_observer, 'ideal')
            GaborData.choice(trial) = GaborData.ideal_answer(trial);
        end

        %% Accuracy & Feedback
        GaborData.accuracy(trial) = GaborData.choice(trial) == GaborData.correct_answer(trial);
        if isempty(GaborData.model_observer)
            if GaborData.accuracy(trial)
                sounds(1, 0.2);
            else
                sounds(0, 0.2);
            end
            pause(0.5); % Pause for 500 ms after feedback before next trial
        end

        trial = trial + 1;
        block_trial = block_trial + 1;
    end
catch ERR
    earlyQuit;
    Screen('CloseAll');
    rethrow(ERR);
end

%% Save final data to folder
Screen('CloseAll');
save(fileName, 'GaborData', 'image_collection');
end