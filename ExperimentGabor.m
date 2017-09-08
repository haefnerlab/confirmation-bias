function GaborData = ExperimentGabor(GaborData, varargin)

directory = fullfile(pwd, '..');
settings = LoadSettings(directory);

datadir = fullfile(directory, 'RawData');
if ~exist(datadir, 'dir'), mkdir(datadir); end

subjectID = getSubjectId(datadir, 'gaborV2');
sessionNo = length(dir(fullfile(datadir, [subjectID '*']))) + 1;
subjectID = [subjectID '-Session' num2str(sessionNo)];

%% Environment and PsychToolBox Initialization

if isempty(GaborData.model_observer)
    
    % Define variables that PTB adds to the 'static workspace' here (avoids an
    % error message...)
    global AGL GL GLU ptb_RootPath ptb_ConfigPath;
    
    cd(fullfile(directory, 'Code')) % Set the current directory
    commandwindow; % Moves the cursor to the commandwindow
    
    if settings.useOpenGL, InitializeMatlabOpenGL; end
    
    % Screen set up
    whichScreen = settings.whichScreen; %allow to choose the display if there's more than one
    xc = settings.screenSize(3)/2; %	Gets the middle of the horizontal axis
    yc = settings.screenSize(4)/2; % Gets the middle of the vertical axis
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
        'fixationSymbol', '+', ...
        'fixationCenter', [xc, yc], ...
        'fixationSymbolSize', [20 20], ...
        'fixationSymbolColors', [180 180 180], ...
        'fixationRadius', 18, ...
        varargin{:});
    
    % Set up keyboard functions
    KbName('UnifyKeyNames');
    goKey = KbName(settings.keyGo);
    exitKey = KbName(settings.keyExit);
    
end

if isequal(GaborData.stair_fn, @Staircase.contrast)
    fileName = fullfile(datadir, [subjectID '-GaborDataContrast.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataContrastQuit.mat']);
elseif isequal(GaborData.stair_fn, @Staircase.ratio)
    fileName = fullfile(datadir, [subjectID '-GaborDataRatio.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataRatioQuit.mat']);
elseif isequal(GaborData.stair_fn, @Staircase.noise)
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

if isempty(GaborData.model_observer)
    % Create 2 textures to display templates.
    right_template = squeeze(bpg.genImages(1, GaborData.stim_size, GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, GaborData.right_category, 2)) * 64.0 + 127.0;
    left_template = squeeze(bpg.genImages(1, GaborData.stim_size, GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, GaborData.left_category, 2)) * 64.0 + 127.0;
    right_tex = Screen('MakeTexture', wPtr, right_template);
    left_tex = Screen('MakeTexture', wPtr, left_template);
    [h, w] = size(right_template);
    
    % Further timing setup
    GaborData.blank_duration = settings.monitorFPS * GaborData.blank_frames;
    
    %% Begin experiment
    EyeTracker.AutoCalibrate(tracker_info);
    
    % Instruction Screen
    textbox = ptbCenteredRect([xc yc/2], settings.screenSize(3:4)/3);
    Screen('TextSize', wPtr, 30); % Set text size to 30
    Screen('FillRect', wPtr, 127);
    DrawFormattedText(wPtr, ...
        ['You will see a series of images flashing very quickly at the center of the screen. ' ...
        'You are required to keep your eyes on the white cross. Then you will be shown two images (examples below). ' ...
        'You will have to decide which image is most consistent with the preceding frames. ' ...
        sprintf('Select the image positioned to the left or right by pressing %s or %s respectively. ', settings.keyLeftName, settings.keyRightName) ...
        'Ask the researcher if you need further clarification. ' ...
        sprintf('Press %s to begin.', settings.keyGoName)], ...
        'centerblock', 'center', white, 100, 0, 0, 1.5, 0, textbox);
    templatey = (textbox(4) + settings.screenSize(4)) / 2;
    Screen('DrawTexture', wPtr, left_tex, [], ptbCenteredRect([xc-w templatey], [w h]));
    Screen('DrawTexture', wPtr, right_tex, [], ptbCenteredRect([xc+w templatey], [w h]));
    Screen('Flip', wPtr); % Function to flip to the next screen image
    if ptbWaitKey([goKey exitKey]) == exitKey
        Screen('CloseAll');
        return;
    end
    Screen('Flip', wPtr); % Function to flip to the next screen image
end

    function earlyQuit
        if GaborData.current_trial > 5
            save(fileNameQuit, 'GaborData');
        end
        ShowCursor();
        Screen('CloseAll');
    end

% Begin Preliminary Trials
seen_block_notification = false;
trial = 1;
block_trial = 1;
block = 1;

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
                if isempty(GaborData.model_observer)
                    sounds(-1, 1.5);
                    DrawFormattedText(wPtr, ...
                        [sprintf('You have completed a %d blocks. ', block) ...
                        'You may take a break if you want! ' ...
                        sprintf('Press %s whenever you are ready again.', settings.keyGoName) ...
                        '\nThe images to be discriminated are displayed again below.'], ...
                        'centerblock', 'center', white, 60, 0, 0, 1.5, 0, textbox);
                    Screen('DrawTexture', wPtr, left_tex, [], ptbCenteredRect([xc-w templatey], [w h]));
                    Screen('DrawTexture', wPtr, right_tex, [], ptbCenteredRect([xc+w templatey], [w h]));
                    Screen('Flip', wPtr);
                    block=block+1;
                    if ptbWaitKey([goKey exitKey]) == exitKey
                        earlyQuit;
                        return;
                    end
                    
                    Screen('Flip', wPtr);
                end
            end
            
            % Start of a block - set params to initial values.
            GaborData.streak(trial) = 0;
            GaborData.reversal_counter(trial) = 0;
            GaborData.contrast(trial) = GaborData.contrast(1);
            GaborData.ratio(trial) = GaborData.ratio(1);
            GaborData.noise(trial) = GaborData.noise(1);
            GaborData.step_size(trial) = GaborData.step_size(1);
            if isfield(GaborData, 'iid')
                GaborData.iid(trial) = GaborData.iid(1);
            end
            block_trial = 1;
        else
            seen_block_notification = false;
            
            GaborData.contrast(trial) = GaborData.contrast(trial-1);
            GaborData.ratio(trial) = GaborData.ratio(trial-1);
            GaborData.noise(trial) = GaborData.noise(trial-1);
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
            
            % Apply the staircase
            GaborData = GaborData.stair_fn(GaborData);
        end
        
        %% Run this trial
        
        % Generate stimulus for this trial.
        [image_array, frame_categories, checksum] = GaborStimulus(GaborData, trial);
        GaborData.frame_categories(trial, :) = frame_categories;
        GaborData.checksum(trial) = checksum;
        
        % Record answer of the ideal observer.
        GaborData.ideal_frame_signals(trial, :) = ...
            bpg.getSignal(image_array - 127, GaborData.left_category, max(GaborData.noise(trial), .04)) - ...
            bpg.getSignal(image_array - 127, GaborData.right_category, max(GaborData.noise(trial), .04));
        GaborData.ideal_answer(trial) = 1 * (sum(GaborData.ideal_frame_signals(trial, :)) > 0);
        
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
        elseif strcmpi(GaborData.model_observer, 'oracle')
            GaborData.choice(trial) = GaborData.correct_answer(trial);
        elseif strcmpi(GaborData.model_observer, 'bernoulli')
            decision_var = dot(GaborData.ideal_frame_signals(trial, :), GaborData.model_pk);
            bernoulli_p = sigmoid(decision_var / GaborData.sigmoid_slope);
            % Choose sign of decision_var with probability related to
            % magnitude of decision_var.
            GaborData.choice(trial) = 1 * (rand < bernoulli_p);
        end
        
        %% Accuracy & Feedback
        GaborData.accuracy(trial) = GaborData.choice(trial) == GaborData.ideal_answer(trial);
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
save(fileName, 'GaborData');
end