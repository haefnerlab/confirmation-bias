function [image_properties, eye_tracker_points, broke_fixation] = trialStimuliGabor(current_trial, image_array, screen, subjectID, Data, automatic, phase, directory, tracker_info)
% trialStimuli displays the animation of several gabor patches in quick
% succession to the subject, or runs through a single trial of the
% experiment.

%screen_frame = Data.screen_frame * .01667; % This is because the refresh rate is 60 Hz
screen_frame = Data.screen_frame * .0083; % This is because the refresh rate is 120 Hz
% This tells me how many ms each image will be on the screen

eye_tracker_points = [];
broke_fixation = false;

%% Make sure to have left/right patch to match the orientations used
%left_default = makeGabor(Data.gabor_sigma_x, Data.gabor_sigma_y, Data.gabor_spatial_frequency, Data.gabor_phase, ...
%    -pi/4, Data.gabor_startpoint, Data.gabor_step, Data.gabor_endpoint);    % Load left template
%right_default = makeGabor(Data.gabor_sigma_x, Data.gabor_sigma_y, Data.gabor_spatial_frequency, Data.gabor_phase, ...
%    pi/4, Data.gabor_startpoint, Data.gabor_step, Data.gabor_endpoint);    % Load right template

res = Data.screen_resolution;

left_default = [[zeros(res,res) zeros(res,res) ones(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) ones(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) ones(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) ones(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) ones(res,res) zeros(res,res) zeros(res,res)]];

right_default = [[zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res)];
    [ones(res,res) ones(res,res) ones(res,res) ones(res,res) ones(res,res)];
    [zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res)];
    [zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res) zeros(res,res)]];

left_patch = (left_default .* 6.0) .* 16.0 + 127.0;   % To give the template a nice large boost for easy viewing
right_patch = (right_default .* 6.0) .* 16.0 + 127.0;   % To give the template a nice large boost for easy viewing

log_regress = zeros(3, Data.number_of_images);

whichScreen = 0; %allow to choose the display if there's more than one
ResolutionScreen = Screen('Resolution', whichScreen); % Gets screen res
ScreenSize = [0 0 ResolutionScreen.width ResolutionScreen.height]; % Sets full screen
xc = ScreenSize(3)/2; %	Gets the middle of the horizontal axis
yc = ScreenSize(4)/2; % Gets the middle of the vertical axis
Screen('Preference', 'SkipSyncTests', 0); % Opens Screen

white = [255 255 255];          % Sets the color to be white
black = [0 0 0];

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
left = KbName('leftArrow');
right = KbName('rightArrow');
up = KbName('upArrow');
down = KbName('downArrow');

try
    
    wPtr = screen;
    % a pointer to refer to the same screen used in the previous trials
    
    left_template = left_default;
    right_template = right_default;
    
    if automatic == 0
        left_template = (left_default .* Data.contrast(current_trial));     % Match the signal:noise to the current image frame
        right_template = (right_default .* Data.contrast(current_trial));
    end
    
    
    [elementX, elementY] = size(left_template);
    flat_left_template = reshape(left_template',[1,elementX*elementY]);     % Transpose to flatten row-by-row
    
    [elementX, elementY] = size(right_template);
    flat_right_template = reshape(right_template',[1,elementX*elementY]);   % Transpose to flatten row-by-row
    
    log_odds = 0;
    
    
    image_texture = zeros(1,Data.number_of_images);
    for i = 1:Data.number_of_images     % create a openGL texture for each image
        image = squeeze(image_array(i,:,:));
        
        % Image Template Difference
        % To get the log odds for each trial to know the best answer
        [elementX, elementY] = size(image);
        flat_image = reshape(image',[1,elementX*elementY]);
        
        image_properties.log_regress(1,i) = dot(flat_image,flat_left_template);
        image_properties.log_regress(2,i) = dot(flat_image,flat_right_template);
        image_properties.log_regress(3,i) = image_properties.log_regress(1,i) - image_properties.log_regress(2,i);
        
        log_odds = log_odds + image_properties.log_regress(3,i);
        
        %{
        if automatic == 0
            image = (image .* 16.0) + 127.0;
            % Scale it into the color range of 0 to 255 if you're going to be showing it to the screen
        end
        %}
        image_texture(i) = Screen('MakeTexture', wPtr, image);
        
    end
    
    image_properties.log_odds = log_odds / Data.number_of_images;
    
    
    
    Screen('FillRect', wPtr, 127.0);        % Make the background gray
    [~, stimOnsetTime] = Screen('Flip', wPtr);
    % Immediately update the display and store the timestamp of the effective update in stimOnsetTime
    
    
    if automatic == 0
        Screen('DrawLine', wPtr, black, xc-res*4, yc-res*4, xc-res*4+15, yc-res*4, 1);     % Show the corners of where the stimulus which is about to appear
        Screen('DrawLine', wPtr, black, xc-res*4, yc-res*4, xc-res*4, yc-res*4+15, 1);
        Screen('DrawLine', wPtr, black, xc-res*4, yc+res*4, xc-res*4+15, yc+res*4, 1);
        Screen('DrawLine', wPtr, black, xc-res*4, yc+res*4, xc-res*4, yc+res*4-15, 1);
        Screen('DrawLine', wPtr, black, xc+res*4, yc-res*4, xc+res*4-15, yc-res*4, 1);
        Screen('DrawLine', wPtr, black, xc+res*4, yc-res*4, xc+res*4, yc-res*4+15, 1);
        Screen('DrawLine', wPtr, black, xc+res*4, yc+res*4, xc+res*4-15, yc+res*4, 1);
        Screen('DrawLine', wPtr, black, xc+res*4, yc+res*4, xc+res*4, yc+res*4-15, 1);

        [is_fixating, tracker_info, eye_tracker_points] = EyeTracker.getFixation(tracker_info, wPtr, 127);
        start_time = eye_tracker_points(1, 1);
        if ~is_fixating
            broke_fixation = true;
            return;
        end
        
        for i = 1:Data.number_of_images
            EyeTracker.drawFixationSymbol(tracker_info, wPtr);
            gaze_point = EyeTracker.getGazePoint(tracker_info, 'pixels');
            if ~EyeTracker.isFixation(tracker_info, gaze_point)
                broke_fixation = true;
                return;
            end
            eye_tracker_points = vertcat(eye_tracker_points, [GetSecs()-start_time gaze_point]);
            Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-900, yc+550, 0);   % Unobtrusive output to screen of the current trial number
            [x_axis, y_axis] = size(image);
            Screen('DrawTexture', wPtr, image_texture(i), [], [xc-x_axis yc-y_axis xc+x_axis yc+y_axis]); %Fill the buffer with the first texture
            Screen('Flip', wPtr);
            [~, stimOnsetTime] = Screen('Flip', wPtr, stimOnsetTime+screen_frame);
            %update the display in screen_frame after the last Flip?
            
            Screen('Close', image_texture(i));% To save on memory and processing time, close the screen afterwards
        end
        
        %if myIsField(Data, 'move_on') %&& Data.current_trial < 26
        show_left_patch = Screen('MakeTexture', wPtr, left_patch);
        Screen('DrawTexture', wPtr, show_left_patch, [], [xc-res*4-200 yc-res*4 xc+res*4-200 yc+res*4]);   % xc, yc indicates the coordinates of the middle of the screen
        show_right_patch = Screen('MakeTexture', wPtr, right_patch);
        Screen('DrawTexture', wPtr, show_right_patch, [], [xc-res*4+200 yc-res*4 xc+res*4+200 yc+res*4]);
        
        
        
        %{
        else
        show_left_patch = Screen('MakeTexture', wPtr, left_patch);
        Screen('DrawTexture', wPtr, show_left_patch, [], [xc-res*5-500 yc-res*5 xc+res*5-500 yc+res*5]);   % xc, yc indicates the coordinates of the middle of the screen
        show_right_patch = Screen('MakeTexture', wPtr, right_patch);
        Screen('DrawTexture', wPtr, show_right_patch, [], [xc-res*5+500 yc-res*5 xc+res*5+500 yc+res*5]);
        end
        %}
        
        
        Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-600, yc+250, 0);   % Unobtrusive output to screen of the current trial number
        onset = Screen('Flip', wPtr);
        
        tstart=tic;
        [~,~,keyCode] = KbCheck;
        while ~keyCode(left) && ~keyCode(right) && toc(tstart)<=1 % wait for press
            [~,~,keyCode] = KbCheck;
            if keyCode(escapeKey)
                
                if ~exist(directory, 'dir')
                    mkdir(directory);
                    
                    fileName = sprintf('%s%s-GaborQuit.mat',directory,subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                else
                    fileName = sprintf('%s%s-GaborQuit.mat',directory,subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                end
                ShowCursor([],whichScreen)
                sca; % closes screen
                return
            end
        end
        
        offset = Screen('Flip', wPtr);
        
        image_properties.reaction = (offset - onset)*1000;  % Records reaction time in ms times a thousand
        
        if keyCode(left)% == KbCheck
            image_properties.choice = 1;        % The subject chose left orientation
        elseif keyCode(right)% == KbCheck
            image_properties.choice = 0;        % The subject chose right orientation
        else
            image_properties.choice = nan;
        end
        
        
    elseif automatic == 1
        % Ideal Observer
        
        % Testing the observer on volume level
        if phase == 0
            threshold = 135;
            
            if Data.contrast(current_trial) >= threshold
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 1;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 0;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            else
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 0;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 1;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 8;
            
            if Data.ratio(current_trial) >= threshold
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 1;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 0;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            else
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 0;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 1;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
        end
        
        image_properties.reaction = 100;
        
    elseif automatic == 2
        % Humanlike Observer with decision noise
        
        % Testing the observer on volume level
        if phase == 0
            threshold = 135;
            sensory_noise = randn * sqrt(15);
            if threshold + sensory_noise < 127
                threshold = 129;
                sensory_noise = 0;
            end
            decision_noise = 0.1;
            
            if Data.contrast(current_trial) >= (threshold + sensory_noise)
                if rand < decision_noise
                    if Data.correct_answer(current_trial) == 0
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                else
                    image_properties.choice = Data.correct_answer(current_trial);
                end
            elseif Data.volume(current_trial) < (threshold + sensory_noise)
                if rand < decision_noise
                    image_properties.choice = Data.correct_answer(current_trial);
                else
                    if Data.correct_answer(current_trial) == 0
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
            
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 8;
            sensory_noise = normrnd(0,0.05);
            if threshold + sensory_noise < Data.number_of_images/2
                threshold = 13;
                sensory_noise = 0;
            end
            decision_noise = 0.001;
            
            if Data.ratio(current_trial) >= (threshold + sensory_noise)
                if rand < decision_noise
                    if Data.correct_answer(current_trial) == 0
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                else
                    image_properties.choice = Data.correct_answer(current_trial);
                end
                
            elseif Data.ratio(current_trial) < (threshold + sensory_noise)
                if rand < decision_noise
                    image_properties.choice = Data.correct_answer(current_trial);
                else
                    if Data.correct_answer(current_trial) == 0
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
        end
        
        image_properties.reaction = 100;
        
    elseif automatic == 3
        % Human Observer
        % Testing the observer on volume level
        if phase == 0
            threshold = 135;
            sensory_noise = 15;
            
            if Data.contrast(current_trial) >= threshold
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 1;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 0;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            else
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 0;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 1;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 8;
            
            if Data.ratio(current_trial) >= threshold
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 1;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 0;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            else
                if sum(Data.order_of_orientations(i,:)) > Data.number_of_images/2
                    image_properties.choice = 0;
                elseif sum(Data.order_of_orientations(i,:)) < Data.number_of_images/2
                    image_properties.choice = 1;
                else
                    if rand < 0.5
                        image_properties.choice = 1;
                    else
                        image_properties.choice = 0;
                    end
                end
            end
        end
        
        image_properties.reaction = 100;
        
    end
    
    
catch ERR
    
    Screen('CloseAll');
    rethrow(ERR);
    
end