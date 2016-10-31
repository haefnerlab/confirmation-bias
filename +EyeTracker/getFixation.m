function [is_holding, tracker_info, coordinates_history] = getFixation(tracker_info, wPtr, background)
%GETFIXATION displays a fixation symbol on a blank background until the
%subject is deemed fixating.
%
% [is_holding, tracker_info, coordinates_history] = GETFIXATION(tracker_info, wPtr, background)
% runs the routine drawing to the screen specified by wPtr. Further
% parameters are controlled by tracker_info. Returns a boolean is_holding,
% an updated tracker_info that has proper slip correction, and the history
% of gaze locations in pixel space with [msec, x, y] in each row.
% 'background' may be a color or a function to draw the background.
start = GetSecs();
coordinates_history = [];
i = 1;
is_holding = false;
while GetSecs() - start < tracker_info.fixationTime / 1000
    if isnumeric(background)
        Screen('FillRect', wPtr, background);
    else
        background();
    end
    EyeTracker.drawFixationSymbol(tracker_info, wPtr);
    now_msecs = (GetSecs() - start) * 1000;
    [gx, gy] = EyeTracker.getGazePoint(tracker_info, 'pixels');
    
    coordinates_history(i,:) = [now_msecs gx gy];
    
    if now_msecs > tracker_info.fixationMinimumHold
        indices_within_minimum_hold = now_msecs - coordinates_history(:,1) < tracker_info.fixationMinimumHold;
        [center, gaze_bbox] = get_bbox(coordinates_history(indices_within_minimum_hold, 2:3));
        bbox_width = gaze_bbox(3) - gaze_bbox(1);
        bbox_height = gaze_bbox(4) - gaze_bbox(2);
        if bbox_width < 2*tracker_info.fixationRadius && bbox_height < 2*tracker_info.fixationRadius
            is_holding = true;
            break;
        end
        if tracker_info.debug
            Screen('FrameRect', wPtr, [255, 0, 0], gaze_bbox);
        end
    end
    if tracker_info.debug
        % draw target area
        Screen('FrameOval', wPtr, [255 255 255], ptbCenteredRect(tracker_info.fixationCenter, 2*tracker_info.fixationRadius*ones(1,2)));
        % draw crosshairs on gaze location
        Screen('DrawLine', wPtr, [255 0 0], gx - 10, gy, gx + 10, gy);
        Screen('DrawLine', wPtr, [255 0 0], gx, gy - 10, gx, gy + 10);
    end

    Screen('Flip', wPtr);
    
    i = i + 1;
end
if is_holding
    tracker_info.fixationCorrection = tracker_info.fixationCenter - center;
end
end

function [center, bbox] = get_bbox(points)
x = points(:, 1);
y = points(:, 2);

center = [mean(x) mean(y)];
bbox = [min(x) min(y) max(x) max(y)];
end