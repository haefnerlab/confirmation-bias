function demoFixation(whichscreen, n_repeats, varargin)

if nargin < 2, n_repeats = 5; end

InitializeMatlabOpenGL;
Screen('Preference', 'SkipSyncTests', 1);
[wPtr, ~] = Screen('OpenWindow', whichscreen, [0 0 0], [], 32);

HideCursor(whichscreen);

tracker_info = EyeTracker.initEyeTracker(whichscreen, varargin{:});
disp(tracker_info);
trajectories = cell(n_repeats, 1);

try
    for n=1:n_repeats
        Screen('FillRect', wPtr, 127.0);
        EyeTracker.drawFixationSymbol(tracker_info, wPtr);
        Screen('Flip', wPtr);
        pause(2);
        [fixating, tracker_info, trajectory] = EyeTracker.getFixation(tracker_info, wPtr, 127);
        disp(fixating);
        trajectories{n} = trajectory(:,2:3);
        disp(trajectory);
        if fixating
            % the task: hold fixation for a few seconds while things happen
            % in the periphery
            start = GetSecs();
            gaze_xy = [];
            while GetSecs() - start < 2
                Screen('FillRect', wPtr, 127);
                rect_corner = floor(rand(1,2).*tracker_info.pixelsPerGazeCoordinate);
                rect = [rect_corner, rect_corner + [50 50]];
                Screen('FillRect', wPtr, randi(256, 1, 3), rect);
                EyeTracker.drawFixationSymbol(tracker_info, wPtr);
                Screen('Flip', wPtr);
                gaze_xy = vertcat(gaze_xy, EyeTracker.getGazePoint(tracker_info));
                if ~EyeTracker.isValidTrial(tracker_info, gaze_xy)
                    Screen('FillRect', wPtr, [255 0 0]);
                    Screen('Flip', wPtr);
                    pause(1);
                    break;
                end
            end
        end
    end
    pause(1);
    
    for i=1:n_repeats
        figure();
        plot(trajectories{i}(:,1), trajectories{i}(:,2));
    end
    
catch ERR
    Screen('CloseAll');
    rethrow(ERR);
end

Screen('CloseAll');
end