function [ success ] = AutoCalibrate(tracker_info)
% AUTOCALIBRATE communicate with the ViewPoint software and run their
% calibration code.
%
% success = AUTOCALIBRATE(tracker_info) uses the following configuration
% options in the tracker_info struct (returning a boolean indicating
% successful calibration):
%
% tracker_info.calibration_n_points - number of points to use; must be in
% {6, 9, 12, 16}
%
% tracker_info.calibration_animate - 'bounce' or 'shrink'
%
% tracker_info.calibration_color - [r g b] ints in [0, 255] specifying the
% stimulus color (default is bright green)

if tracker_info.debugUseMouse
    success = true;
    return;
end

assert(any(tracker_info.calibration_n_points == [6 9 12 16]), 'calibration_n_points must be in {6, 9, 12, 16}');

% Psychtoolbox can't be open while calibration runs.
Screen('CloseAll');

% Send configuration commands.
vpx_SendCommandString(sprintf('calibration_Points %d', tracker_info.calibration_n_points));
vpx_SendCommandString('calibration_PointLocationMethod Automatic');
vpx_SendCommandString(sprintf('calibration_StimulusType %s', tracker_info.calibration_animate));
r = tracker_info.calibration_color(1);
g = tracker_info.calibration_color(2);
b = tracker_info.calibration_color(3);
vpx_SendCommandString(sprintf('calibration_StimulusColor %d %d %d', r, g, b));

% Run calibration.
vpx_SendCommandString('calibrationStart');
vpx_CalibrationInProgress = 6; % See documentation for vpx_GetStatus.
pause(0.1);
while vpx_GetStatus(vpx_CalibrationInProgress), end
success = true;
% TODO - validate result with vpx_getcalibrationeventrecord ??

Screen('OpenWindow', tracker_info.whichscreen);
HideCursor(tracker_info.whichscreen);
end