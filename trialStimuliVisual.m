function image_properties = trialStimuliVisual(screen, subjectID, Data, automatic, directory, settings)
% trialStimuli displays the animation of several gabor patches in quick
% succession to the subject, or runs through a single trial of the
% experiment.


directory = [directory 'RawData'];    % Directory to save the data and files to


screen_frame = Data.screen_frame * .01667; % This is because the refresh rate is 60 Hz
% This tells me how many ms each image will be on the screen

whichScreen = settings.whichScreen; %allow to choose the display if there's more than one
ResolutionScreen = Screen('Resolution', whichScreen); % Gets screen resolution
ScreenSize = [0 0 ResolutionScreen.width ResolutionScreen.height]; % Sets full screen
xc = ScreenSize(3)/2; %	Gets the middle of the horizontal axis
yc = ScreenSize(4)/2; % Gets the middle of the vertical axis

white = [255 255 255];          % Sets the color to be white
black = [0 0 0];

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
exitKey = KbName(settings.keyExit);
leftKey = KbName(settings.keyLeft);
rightKey = KbName(settings.keyRight);

try
    
    wPtr = screen;
    % a pointer to refer to the same screen used in the previous trials
    
    x_axis_pixels = Data.image_length_x;
    y_axis_pixels = Data.image_length_y;
    
    resolution_x = Data.screen_resolution*(x_axis_pixels/2);
    resolution_y = Data.screen_resolution*(y_axis_pixels/2);
    
    stimulus_separation = 250;  % How many pixels to keep between the two stimuli?
    
    
    if myIsField(Data,'move_on')
        screen_frame = .2;   % The bar will only appear for 200 ms in the preliminary phase since it's a bar detection task
    end
    
    image_texture = zeros(2,Data.number_of_images);
    for i = 1:Data.number_of_images     % create a openGL texture for each image
        left_color_value = Data.order_of_flashes(Data.current_trial,1,i);
        right_color_value = Data.order_of_flashes(Data.current_trial,2,i);
        
        left_image = left_color_value*ones(Data.image_length_y,Data.image_length_x);
        right_image = right_color_value*ones(Data.image_length_y,Data.image_length_x);
        
        image_texture(1,i) = Screen('MakeTexture', wPtr, left_image);
        image_texture(2,i) = Screen('MakeTexture', wPtr, right_image);
    end
    
    
    
    Screen('FillRect', wPtr, 127.5);        % Make the background gray
    [~, stimOnsetTime] = Screen('Flip', wPtr);
    % Immediately update the display and store the timestamp of the effective update in stimOnsetTime
    
    
    if automatic == 0
        
        Screen('DrawDots', wPtr, [xc, yc], 50, black, [0 0], 1);
        Screen('DrawDots', wPtr, [xc, yc], 15, white, [0 0], 1);   % draw a bull's eye below where the stimulus will appear as a fixation point
        Screen('DrawLine', wPtr, white, xc-25, yc, xc+25, yc, 1);
        Screen('DrawLine', wPtr, white, xc, yc-25, xc, yc+25, 1);   % draw a cross on the bull's eye as a fixation point
        
        Screen('DrawLine', wPtr, black, xc-resolution_x-stimulus_separation, yc-resolution_y, xc-resolution_x-stimulus_separation+15, yc-resolution_y, 1);     % Show the corners of where the stimulus which is about to appear
        Screen('DrawLine', wPtr, black, xc-resolution_x-stimulus_separation, yc-resolution_y, xc-resolution_x-stimulus_separation, yc-resolution_y+15, 1);
        Screen('DrawLine', wPtr, black, xc-resolution_x-stimulus_separation, yc+resolution_y, xc-resolution_x-stimulus_separation+15, yc+resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc-resolution_x-stimulus_separation, yc+resolution_y, xc-resolution_x-stimulus_separation, yc+resolution_y-15, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x-stimulus_separation, yc-resolution_y, xc+resolution_x-stimulus_separation-15, yc-resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x-stimulus_separation, yc-resolution_y, xc+resolution_x-stimulus_separation, yc-resolution_y+15, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x-stimulus_separation, yc+resolution_y, xc+resolution_x-stimulus_separation-15, yc+resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x-stimulus_separation, yc+resolution_y, xc+resolution_x-stimulus_separation, yc+resolution_y-15, 1);
        
        Screen('DrawLine', wPtr, black, xc-resolution_x+stimulus_separation, yc-resolution_y, xc-resolution_x+stimulus_separation+15, yc-resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc-resolution_x+stimulus_separation, yc-resolution_y, xc-resolution_x+stimulus_separation, yc-resolution_y+15, 1);
        Screen('DrawLine', wPtr, black, xc-resolution_x+stimulus_separation, yc+resolution_y, xc-resolution_x+stimulus_separation+15, yc+resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc-resolution_x+stimulus_separation, yc+resolution_y, xc-resolution_x+stimulus_separation, yc+resolution_y-15, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x+stimulus_separation, yc-resolution_y, xc+resolution_x+stimulus_separation-15, yc-resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x+stimulus_separation, yc-resolution_y, xc+resolution_x+stimulus_separation, yc-resolution_y+15, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x+stimulus_separation, yc+resolution_y, xc+resolution_x+stimulus_separation-15, yc+resolution_y, 1);
        Screen('DrawLine', wPtr, black, xc+resolution_x+stimulus_separation, yc+resolution_y, xc+resolution_x+stimulus_separation, yc+resolution_y-15, 1);
        Screen('Flip', wPtr);
        pause(0.4)
        Screen('Flip', wPtr);
        
        for i = 1:Data.number_of_images
            
            Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-900, yc+550, 0);   % Unobtrusive output to screen of the current trial number
            
            Screen('DrawTexture', wPtr, image_texture(1,i), [], [xc-resolution_x-stimulus_separation yc-resolution_y xc+resolution_x-stimulus_separation yc+resolution_y]); %Draw the bars
            Screen('DrawTexture', wPtr, image_texture(2,i), [], [xc-resolution_x+stimulus_separation yc-resolution_y xc+resolution_x+stimulus_separation yc+resolution_y]);
            
            Screen('DrawDots', wPtr, [xc, yc], 50, black, [0 0], 1);
            Screen('DrawDots', wPtr, [xc, yc], 15, white, [0 0], 1);   % draw a bull's eye below/between where the stimulus will appear as a fixation point
            Screen('DrawLine', wPtr, white, xc-25, yc, xc+25, yc, 1);
            Screen('DrawLine', wPtr, white, xc, yc-25, xc, yc+25, 1);   % draw a cross on the bull's eye as a fixation point
            [~, stimOnsetTime] = Screen('Flip', wPtr, stimOnsetTime+screen_frame);
            %update the display in screen_frame after the last Flip?
            Screen('Close', image_texture(1,i)); % To save on memory and processing time, close the screen afterwards
            Screen('Close', image_texture(2,i));
        end
        
        if Data.current_trial < 26
            Screen('DrawText', wPtr, 'From which side did you see more bright flashes?', xc-500, yc-150, black);
            
            Screen('DrawLine', wPtr, black, xc-resolution_x-200, yc, xc-resolution_x-400, yc, 10);         % Draw the line of the arrow
            Screen('DrawLine', wPtr, black, xc+resolution_x+200, yc, xc+resolution_x+400, yc, 10);
            Screen('DrawLine', wPtr, black, xc-resolution_x-400, yc+15, xc-resolution_x-400, yc-15, 2);    % Draw the head of the arrow
            Screen('DrawLine', wPtr, black, xc+resolution_x+400, yc+15, xc+resolution_x+400, yc-15, 2);
            Screen('DrawLine', wPtr, black, xc-resolution_x-400, yc-15, xc-resolution_x-400-50, yc, 2);
            Screen('DrawLine', wPtr, black, xc+resolution_x+400, yc-15, xc+resolution_x+400+50, yc, 2);
            Screen('DrawLine', wPtr, black, xc-resolution_x-400, yc+15, xc-resolution_x-400-50, yc, 2);
            Screen('DrawLine', wPtr, black, xc+resolution_x+400, yc+15, xc+resolution_x+400+50, yc, 2);
        end
        Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-900, yc+550, 0);   % Unobtrusive output to screen of the current trial number
        onset = Screen('Flip', wPtr);
        
        [~,~,keyCode] = KbCheck;
        while ~keyCode(leftKey) && ~keyCode(rightKey) % wait for press
            [~,~,keyCode] = KbCheck;
            if keyCode(exitKey)
                
                if ~exist(directory, 'dir')
                    mkdir(directory);
                    
                    fileName = sprintf('%s%s-VisualQuit.mat',directory,subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                else
                    fileName = sprintf('%s%s-VisualQuit.mat',directory,subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                end
                sca; % closes screen
                return
            end
        end
        
        offset = Screen('Flip', wPtr);
        
        image_properties.reaction = (offset - onset)*1000;  % Records reaction time in ms times a thousand
        
        if keyCode(leftKey) == KbCheck
            image_properties.choice = 1;        % The subject chose the left image
        elseif keyCode(rightKey) == KbCheck
            image_properties.choice = 0;        % The subject chose the right image
        end
        
    elseif automatic == 1
        answer = Data.correct_answer(Data.current_trial);
        if answer == 1
            image_properties.choice = 1;        % The subject chose the left image
        elseif answer == 0
            image_properties.choice = 0;        % The subject chose the right image
        else
            if rand < 0.5
                image_properties.choice = 1;        % The subject chose the left image
            else
                image_properties.choice = 0;        % The subject chose the right image
            end
        end
        image_properties.reaction = 100;
    end
    
    
catch ERR
    
    Screen('CloseAll');
    rethrow(ERR);
    
end