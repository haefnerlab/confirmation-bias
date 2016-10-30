function stimulus_properties = trialStimuliAuditory(screen, subjectID, Data, current_trial, automatic, phase, add_noise, directory, settings)
% trialStimuli will play the sound several clicks into both ears,
% headphones are required

% Ex. trialStimuliAuditory(0, 'subject', Test_Data, 50, 50, 0, pwd, LoadSettings(pwd))


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

try    % If there is ever an error during an experiment, the PsychToolBox screen will automatically exit in the catch block
    
    wPtr = screen;
    % a pointer to refer to the same screen used in the previous trials
    % always 0
    if automatic ~=5
        Screen('FillRect', wPtr, 127.0);        % Make the background gray
        [~, stimOnsetTime] = Screen('Flip', wPtr);
        % Immediately update the display and store the timestamp of the effective update in stimOnsetTime
    end
    if automatic == 0
        
        % Sum 2, 4, 6, 8, and 16 kHz pure tones together
        tone = tone_generator(6000, 3, 1, 2000, pi/3, 1, @(N)(hanning(N).^2)) + ...
            tone_generator(6000, 3, 1, 4000, pi/3, 1, @(N)(hanning(N).^2)) + ...
            tone_generator(6000, 3, 1, 6000, pi/3, 1, @(N)(hanning(N).^2)) + ...
            tone_generator(6000, 3, 1, 8000, pi/3, 1, @(N)(hanning(N).^2)) + ...
            tone_generator(6000, 3, 1, 16000, pi/3, 1, @(N)(hanning(N).^2));
        cosine = cos(linspace(0,2*pi,19));
        click_sound = conv(tone, cosine);  % Convolve the sum of tones with cosine
        
        maxElement = max(click_sound);    % Rescale the click into a -1 to 1 scaling (the sound vector has to be within that range)
        minElement = min(click_sound);
        click_sound = ((click_sound - minElement)./(maxElement-minElement) - 0.5 ) * 2; % Click sound created!
        
        
        [number_of_left_clicks, number_of_right_clicks, clicks] = makePoissonClicks(Data.number_of_frames, Data.average_clicks(1,current_trial), ...
            Data.average_clicks(2,current_trial), Data.volume(current_trial), Data.volume(current_trial));
        % Generate a 120-frame vector of 1s and 0s to denote a click or lack of click per click-frame
        y = zeros(2,Data.bins);
        bin1 = Data.sampling_rate*0.1;
        frames = 6;
        t = zeros(2,bin1);
        r = zeros(2,frames);
        for i=1:size(r,2)
            if mod(i,2) == 1
                r(:,i)=1;
            end
        end
        for i = 1:length(clicks)
            bin_interval = Data.bins/Data.number_of_frames;   % Determine number of vector elements each bin will have
            if clicks(1,i) == Data.volume(current_trial)    % Embed a click in the middle of the bin
                y(1, ceil(bin_interval*(i-1) + (bin_interval/2)-(length(click_sound)/2)):ceil(bin_interval*(i-1) + (bin_interval/2)+(length(click_sound)/2))-1) = click_sound*Data.volume(current_trial);
            end
            if clicks(2,i) == Data.volume(current_trial)    % Embed a click in the middle of the bin
                y(2, ceil(bin_interval*(i-1) + (bin_interval/2)-(length(click_sound)/2)):ceil(bin_interval*(i-1) + (bin_interval/2)+(length(click_sound)/2))-1) = click_sound*Data.volume(current_trial);
            end
        end
        if add_noise == 1
            y = y + 0.0015*randn(2, length(y));  % Add a little white noise
        end
        
        for i=1:size(r, 2)
            bin_interval= bin1/frames;   % Determine number of vector elements each bin will have
            if r(1,i) == 1    % Embed a click in the middle of the bin
                t(1, ceil(bin_interval*(i-1) + (bin_interval/2)-(length(click_sound)/2)):ceil(bin_interval*(i-1) + (bin_interval/2)+(length(click_sound)/2))-1) = click_sound*Data.volume(current_trial);
            end
            if r(2,i) == 1    % Embed a click in the middle of the bin
                t(2, ceil(bin_interval*(i-1) + (bin_interval/2)-(length(click_sound)/2)):ceil(bin_interval*(i-1) + (bin_interval/2)+(length(click_sound)/2))-1) = click_sound*Data.volume(current_trial);
            end
            
        end
        sound(t,Data.sampling_rate)
        pause(.1);
        sound(y, Data.sampling_rate)         % Play sound
        pause(Data.stimulus_duration)
        stimulus_properties.number_of_left_clicks = number_of_left_clicks;   % Save number of clicks playing per ear
        stimulus_properties.number_of_right_clicks = number_of_right_clicks;
        stimulus_properties.clicks = clicks;   % Save the binary vector denoting the prescence or lack of clickes
        
        %       if myIsField(Data, 'step_size') %&& current_trial < 26     % For the first 25 trials, put instructions telling the subject to choose left or right
        
        Screen('DrawText', wPtr, 'From which direction did you hear more clicks?', xc-500, yc-150, black);
        
        Screen('DrawLine', wPtr, black, xc-200, yc, xc-400, yc, 10);         % Draw the line of the arrow
        Screen('DrawLine', wPtr, black, xc+200, yc, xc+400, yc, 10);
        Screen('DrawLine', wPtr, black, xc-400, yc+15, xc-400, yc-15, 2);    % Draw the head of the arrow
        Screen('DrawLine', wPtr, black, xc+400, yc+15, xc+400, yc-15, 2);
        Screen('DrawLine', wPtr, black, xc-400, yc-15, xc-400-50, yc, 2);
        Screen('DrawLine', wPtr, black, xc+400, yc-15, xc+400+50, yc, 2);
        Screen('DrawLine', wPtr, black, xc-400, yc+15, xc-400-50, yc, 2);
        Screen('DrawLine', wPtr, black, xc+400, yc+15, xc+400+50, yc, 2);
        %    end
        %% Need to edit this for second preliminary phase to show proper number
        Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-550, yc+325, 0);   % Unobtrusive output to screen of the current trial number
        onset = Screen('Flip', wPtr);
        
        
        %         [KeyIsDown,~,keyCode] = KbCheck;
        %         while ~KeyIsDown
        %             [KeyIsDown, ~, keyCode]=KbCheck;
        % set up the timer
        %tt = timer ;
        %set(0,'DefaultFigureVisible','off');
        %tt.timerfcn = 'uiresume' ;
        %tt.startdelay = 2 ;
        %start(tt);
        tstart=tic;
        [~, keyCode, ~] = KbPressWait([], 2);
        while ~keyCode(leftKey) && ~keyCode(rightKey) && toc(tstart)<=1 %strcmp(get(tt,'Running'),'on') % wait for press
            %             [~,~,keyCode] = KbCheck;
            [~, keyCode, ~] = KbPressWait([], 2);
            
            if keyCode(exitKey)
                
                if ~exist([directory 'RawData/'], 'dir')
                    mkdir([directory 'RawData/']);
                    
                    fileName = sprintf('%s%s-AuditoryQuit.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                else
                    fileName = sprintf('%s%s-AuditoryQuit.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
                    save(fileName, 'Data'); % save the data
                end
                sca; % closes screen
                return
            end
        end
        
        % clean up the timer ...
        %stop(tt) ;
        %delete(tt) ;
        
        %         end
        
        %         while ~keyCode(left) && ~keyCode(right) % wait for press
        %             [~,~,keyCode] = KbCheck;
        %             if keyCode(escapeKey)
        %
        %                 if ~exist([directory 'RawData/'], 'dir')
        %                     mkdir([directory 'RawData/']);
        %
        %                     fileName = sprintf('%s%s-AuditoryQuit.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
        %                     save(fileName, 'Data'); % save the data
        %                 else
        %                     fileName = sprintf('%s%s-AuditoryQuit.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save as a csv
        %                     save(fileName, 'Data'); % save the data
        %                 end
        %                 sca; % closes screen
        %                 return
        %             end
        %         end
        
        offset = Screen('Flip', wPtr);
        
        stimulus_properties.reaction = (offset - onset)*1000;  % Records reaction time in ms times a thousand
        
        if keyCode(leftKey)% == KbCheck
            stimulus_properties.choice = 1;        % The subject chose left orientation
        elseif keyCode(rightKey)% == KbCheck
            stimulus_properties.choice = 0;        % The subject chose right orientation
        else
            stimulus_properties.choice = nan;
        end
        
        if number_of_left_clicks > number_of_right_clicks    % The left side has more clicks
            stimulus_properties.answer = 1;
        elseif number_of_left_clicks < number_of_right_clicks    % The right side has more clicks
            stimulus_properties.answer = 0;
        else
            if rand < 0.5    % Both sides has the same number of clicks
                stimulus_properties.answer = 1;
            else
                stimulus_properties.answer = 0;
            end
        end
        
    elseif automatic == 1
        % Ideal Observer
        
        [number_of_left_clicks, number_of_right_clicks, clicks] = makePoissonClicks(Data.number_of_frames, Data.average_clicks(1,current_trial), ...
            Data.average_clicks(2,current_trial), Data.volume(current_trial), Data.volume(current_trial));
        
        % Testing the observer on volume level
        if phase == 0
            threshold = 0.01;
            
            if Data.volume(current_trial) >= threshold
                if number_of_left_clicks > number_of_right_clicks
                    stimulus_properties.choice = 1;
                elseif number_of_left_clicks < number_of_right_clicks
                    stimulus_properties.choice = 0;
                else
                    if rand < 0.5
                        stimulus_properties.choice = 0;
                    else
                        stimulus_properties.choice = 1;
                    end
                end
                
            elseif Data.volume(current_trial) < threshold
                
                if rand < 0.5
                    stimulus_properties.choice = 0;
                else
                    stimulus_properties.choice = 1;
                end
                
            end
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 15;
            
            if Data.ratio(current_trial) >= threshold
                if number_of_left_clicks > number_of_right_clicks
                    stimulus_properties.choice = 1;
                elseif number_of_left_clicks < number_of_right_clicks
                    stimulus_properties.choice = 0;
                else
                    if rand < 0.5
                        stimulus_properties.choice = 0;
                    else
                        stimulus_properties.choice = 1;
                    end
                end
                
            elseif Data.ratio(current_trial) < threshold
                if number_of_left_clicks > number_of_right_clicks
                    stimulus_properties.choice = 0;
                elseif number_of_left_clicks < number_of_right_clicks
                    stimulus_properties.choice = 1;
                else
                    if rand < 0.5
                        stimulus_properties.choice = 0;
                    else
                        stimulus_properties.choice = 1;
                    end
                end
            end
        end
        
        stimulus_properties.reaction = 100;
        stimulus_properties.number_of_left_clicks = number_of_left_clicks;
        stimulus_properties.number_of_right_clicks = number_of_right_clicks;
        stimulus_properties.clicks = clicks;
        
    elseif automatic == 2
        % Humanlike Observer with decision noise
        
        %% Is correct answer based on underlying click rate or actual number of clicks?
        
        [number_of_left_clicks, number_of_right_clicks, clicks] = makePoissonClicks(Data.number_of_frames, Data.average_clicks(1,current_trial), ...
            Data.average_clicks(2,current_trial), Data.volume(current_trial), Data.volume(current_trial));
        
        % Testing the observer on volume level
        if phase == 0
            threshold = 0.01;
            sensory_noise = randn * sqrt(0.5);
            if threshold + sensory_noise < 0
                threshold = 0.001;
                sensory_noise = 0;
            end
            decision_noise = 0.1;
            
            if Data.volume(current_trial) >= (threshold + sensory_noise)
                if rand < decision_noise
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                else
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                end
            elseif Data.volume(current_trial) < (threshold + sensory_noise)
                if rand < decision_noise
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                else
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                end
            end
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 15;
            sensory_noise = normrnd(0,0.05);
            if threshold + sensory_noise < 12
                threshold = 13;
                sensory_noise = 0;
            end
            decision_noise = 0.001;
            
            if Data.ratio(current_trial) >= (threshold + sensory_noise)
                if rand < decision_noise
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                else
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                end
                
            elseif Data.ratio(current_trial) < (threshold + sensory_noise)
                if rand < decision_noise
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                else
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                end
            end
        end
        
        
        
        stimulus_properties.reaction = 100;
        stimulus_properties.number_of_left_clicks = number_of_left_clicks;
        stimulus_properties.number_of_right_clicks = number_of_right_clicks;
        stimulus_properties.clicks = clicks;
        
    elseif automatic == 3
        % Human Observer
        
        %% Is correct answer based on underlying click rate or actual number of clicks?
        
        [number_of_left_clicks, number_of_right_clicks, clicks] = makePoissonClicks(Data.number_of_frames, Data.average_clicks(1,current_trial), ...
            Data.average_clicks(2,current_trial), Data.volume(current_trial), Data.volume(current_trial));
        
        % Testing the observer on volume level
        if phase == 0
            
            % Generate a 120-frame vector of 1s and 0s to denote a click or lack of click per click-frame
            threshold = 0.01;
            sensory_noise = 0.1;
            clicks_volume = Data.volume(current_trial).*clicks + normrnd(0,sensory_noise,[2, length(clicks)]);
            decision_variable = zeros(2);
            for i = 1:Data.number_of_frames
                decision_variable(1) = decision_variable(1) + ((clicks_volume(1,i) - threshold)>=0);
                decision_variable(2) = decision_variable(2) + ((clicks_volume(2,i) - threshold)>=0);
            end
            
            if decision_variable(1)> decision_variable(2)
                stimulus_properties.choice = 1;
            elseif decision_variable(1)< decision_variable(2)
                stimulus_properties.choice = 0;
            else
                if rand < 0.5
                    stimulus_properties.choice = 1;
                else
                    stimulus_properties.choice = 0;
                end
            end
            
            stimulus_properties.reaction = 100;
            stimulus_properties.number_of_left_clicks = number_of_left_clicks;
            stimulus_properties.number_of_right_clicks = number_of_right_clicks;
            stimulus_properties.clicks = clicks;
            
            
            
            
            % Testing the observer on ratio level
        elseif phase == 1
            threshold = 15;
            sensory_noise = normrnd(0,0.5);
            if threshold + sensory_noise < 12
                threshold = 13;
                sensory_noise = 0;
            end
            decision_noise = 0.001;
            
            if Data.ratio(current_trial) >= (threshold + sensory_noise)
                if rand < decision_noise
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                else
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                end
                
            elseif Data.ratio(current_trial) < (threshold + sensory_noise)
                if rand < decision_noise
                    stimulus_properties.choice = Data.correct_answer(current_trial);
                else
                    if Data.correct_answer(current_trial) == 0
                        stimulus_properties.choice = 1;
                    else
                        stimulus_properties.choice = 0;
                    end
                end
            end
            
        end
        
        
        stimulus_properties.reaction = 100;
        stimulus_properties.number_of_left_clicks = number_of_left_clicks;
        stimulus_properties.number_of_right_clicks = number_of_right_clicks;
        stimulus_properties.clicks = clicks;
        set(0,'DefaultFigureVisible','on');
        
        
    elseif automatic == 5
        % Ideal Observer
        
        [number_of_left_clicks, number_of_right_clicks, clicks,l,r] = makePoissonClicksIdeal(Data.number_of_frames, Data.average_clicks(1,current_trial), ...
            Data.average_clicks(2,current_trial), Data.volume(current_trial), Data.volume(current_trial));
        
        % Testing the observer on volume level
          
        if number_of_left_clicks > number_of_right_clicks
            stimulus_properties.choice = 1;
        elseif number_of_left_clicks < number_of_right_clicks
            stimulus_properties.choice = 0;
        else
            if rand < 0.5
                stimulus_properties.choice = 0;
            else
                stimulus_properties.choice = 1;
            end
        end
        
        stimulus_properties.reaction = 100;
        stimulus_properties.number_of_left_clicks = number_of_left_clicks;
        stimulus_properties.number_of_right_clicks = number_of_right_clicks;
        stimulus_properties.clicks = clicks;
        stimulus_properties.l = l;
        stimulus_properties.r = r;
        
    end
catch ERR
    
    Screen('CloseAll');
    rethrow(ERR);
    
end