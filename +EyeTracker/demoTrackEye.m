function demoTrackEye(whichscreen, length_tail, varargin)

InitializeMatlabOpenGL;

[wPtr, ~] = Screen('OpenWindow', whichscreen, [0 0 0], [], 32);

HideCursor(whichscreen);

tracker_info = EyeTracker.initEyeTracker(whichscreen, varargin{:});

history = zeros(length_tail, 2);
i = 1;

tstart = GetSecs();

try
    while GetSecs() - tstart < 10
        Screen('FillRect', wPtr, 127);
        gxgy = EyeTracker.getGazePoint(tracker_info);
        history(i, :) = gxgy;
        for j=1:length_tail
            if j == i, continue; end
            k = mod(j, length_tail) + 1;
            Screen('DrawLine', wPtr, [1 0 0], history(j, 1), history(j, 2), history(k, 1), history(k, 2));
        end
        Screen('Flip', wPtr);
        i = mod(i, length_tail) + 1;
    end
catch ERR
    Screen('CloseAll');
    rethrow(ERR);
end

Screen('CloseAll');
end