function tracker_info = initEyeTracker(whichscreen, varargin)
% Connect to eye tracker.
vpx_Initialize();

% Return struct with tracker's configuration options.
resolution = Screen('Resolution', whichscreen);

% This struct defines all configurable options for the eye tracker code.
tracker_info = struct(...
    'whichscreen', whichscreen, ...
    'debug', false, ...
    'debugUseMouse', false, ...
    'pixelsPerGazeCoordinate', [resolution.width, resolution.height], ... % X, Y screen pixels per 'gaze unit'
    'fixationSymbol', 'r', ... % 'r' for rect, 'c' for circle, or '+' for plus
    'fixationSymbolSize', [10, 10], ... % pixel size of fixation symbol, independent of the 'Rect' below
    'fixationTime', 1000, ... % ms. Max time allowed in getFixation()
    'fixationMinimumHold', 500, ... % Time required within fixationRect to consider it held.
    'fixationCorrection', [0 0], ... % Add this to [gx, gy] to get corrected position (recalibrated during getFixation)
    'fixationCenter', [resolution.width/2, resolution.height/2], ...
    'fixationRect', [30 30]);

% Parse any extra (..., 'key', value, ...) pairs passed in through
% varargin.
for val_idx=2:2:length(varargin)
    key = varargin{val_idx-1};
    if ~ischar(key)
        warning('invalid input to initEyeTracker. After whichscreen, all arguments should be (..., ''key'', value, ...)');
    elseif ~isfield(tracker_info, key)
        warning('unrecognized tracker_info field: ''%s''', key);
    else
        tracker_info.(key) = varargin{val_idx};
    end
end
end