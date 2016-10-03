function [] = Bar_Demo(flash_gap)

% A function to create a bar in the center of the screen and the bar will
% keep flashing until the viewer presses the spacebar.
% flash_gap determines how long it will be between flashes
% Ex. flash_gap = 1 means flash, no flash, flash, no flash
%     flash_gap = 2 means flash, no flash, no flash, flash, no flash, no flash


commandwindow; % Moves the cursor to the commandwindow

InitializeMatlabOpenGL

whichScreen = 0; %allow to choose the display if there's more than one
ResolutionScreen = Screen('Resolution', whichScreen); % Gets screen resolution
ScreenSize = [0 0 ResolutionScreen.width ResolutionScreen.height]; % Sets full screen
xc = ScreenSize(3)/2; %	Gets the middle of the horizontal axis
yc = ScreenSize(4)/2; % Gets the middle of the vertical axis
Screen('Preference', 'SkipSyncTests', 0); % Opens Screen

white = [255 255 255];          % Sets the color to be white
black = [0 0 0];

[wPtr, ~] = Screen('OpenWindow', whichScreen, black, [], 32);

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
left = KbName('leftArrow');
right = KbName('rightArrow');
up = KbName('upArrow');
down = KbName('downArrow');



try
    
    Screen('FillRect', wPtr, 127.0);        % Make the background gray
    [~, stimOnsetTime] = Screen('Flip', wPtr);
    
    color_values = [255 127*ones(1,flash_gap)];
    image_texture = zeros(length(color_values));
    
    [~,~,keyCode] = KbCheck;
    while ~keyCode(spaceKey) % wait for press
        [~,~,keyCode] = KbCheck;
        
        for i = 1:length(color_values)
            image = color_values(i);
            image_texture(i) = Screen('MakeTexture', wPtr, image);
        end
        for i = 1:length(color_values)
            Screen('DrawTexture', wPtr, image_texture(i), [], [xc-25 yc-250 xc+25 yc+250]); % Draw the bars
            [~, stimOnsetTime] = Screen('Flip', wPtr, stimOnsetTime+0.01667);
            %update the display in screen_frame after the last Flip?
            Screen('Close', image_texture(i));
        end
        
        if keyCode(escapeKey)
            sca; % closes screen
            return
        end
    end
    sca;

catch ERR
    
    Screen('CloseAll');
    rethrow(ERR);
    
end

end