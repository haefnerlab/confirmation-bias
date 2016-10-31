function tracker_info = initEyeTracker(whichscreen, varargin)
% Return struct with tracker's configuration options.
resolution = Screen('Resolution', whichscreen);

% This struct defines all configurable options for the eye tracker code.
tracker_info = struct(...
    ... % PTB info and debug flags
    'whichscreen', whichscreen, ...
    'debug', false, ... % when true, renders additional info onto the screen.
    'debugUseMouse', false, ... % when true, substitutes the mouse coordinate for the eye tracker.
    'pixelsPerGazeCoordinate', [resolution.width, resolution.height], ... % X, Y screen pixels per 'gaze unit'
    ... % parameters for getFixation()
    'fixationSymbol', 'r', ... % 'r' for rect, 'c' for circle, 'b' for bullseye, or '+' for plus
    'fixationSymbolSize', [10, 10], ... % pixel size of fixation symbol, independent of the 'Rect' below
    'fixationSymbolColors', [255 255 255; 0 0 0], ... % primary/secondary color of fixation symbol
    'fixationTime', 1000, ... % ms. Max time allowed in getFixation()
    'fixationMinimumHold', 500, ... % Time required within fixation area to consider it held.
    ... % parameters for isFixation()
    'fixationCorrection', [0 0], ... % Add this to [gx, gy] to get corrected position (this is set automatically during getFixation)
    'fixationCenter', [resolution.width/2, resolution.height/2], ...
    'fixationRadius', 10, ... % true size for fixation requirement (separate from the symbol size above)
    ... % parameters for smoothing
    'smoothing_n_points', 4, ... % how far in the past to average together
    ... % parameters for calibration
    'calibration_n_points', 12, ... % number of points used for calibration
    'calibration_animate', 'Shrink', ... % 'Shrink' or 'Bounce'
    'calibration_color', [100 255 100]);

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

if ~tracker_info.debugUseMouse
    % Connect to eye tracker and do one-time configuration commands
    vpx_Initialize();
    vpx_SendCommandString(sprintf('smoothingPoints %d', tracker_info.smoothing_n_points));
end
end