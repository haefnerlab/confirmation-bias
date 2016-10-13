function varargout = getGazePoint(tracker_info, units)
%GETGAZEPOINT get xy location of gaze from eye tracker.
%
% xy = GETGAZEPOINT(tracker_info) uses screen size info in tracker_info to
% return the x and y coordinate of the gaze in pixels.
%
% [x, y] = GETGAZEPOINT(tracker_info) is identical to the above but
% separates out x and y in the returned values, for convenience.
%
% [x, y] = GETGAZEPOINT(tracker_info, 'gaze') gets the coordinates in 'gaze
% space' (i.e. in [0,1] along both x and y) rather than in pixel space.
%
% For debugging, if tracker_info.debugUseMouse is true, this function
% returns the location of the mouse instead of the location of the eyes.
if nargin < 2, units = 'pixels'; end
if tracker_info.debugUseMouse
    [gx, gy] = GetMouse(tracker_info.whichscreen);
    if strcmpi(units, 'gaze')
        gx = gx / tracker_info.pixelsPerGazeCoordinate(1);
        gy = gy / tracker_info.pixelsPerGazeCoordinate(2);
    end
else
    [gx, gy] = vpx_GetGazePoint();
    if strcmpi(units, 'pixels')
        gx = gx * tracker_info.pixelsPerGazeCoordinate(1);
        gy = gy * tracker_info.pixelsPerGazeCoordinate(2);
    end
end
if nargout == 1
    varargout{1} = [gx, gy];
elseif nargout == 2
    varargout{1} = gx;
    varargout{2} = gy;
end
end