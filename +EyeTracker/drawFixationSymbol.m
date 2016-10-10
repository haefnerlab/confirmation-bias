function drawFixationSymbol(tracker_info, wPtr, color)
if nargin < 2, color = [255 255 255]; end
xc = tracker_info.fixationCenter(1);
yc = tracker_info.fixationCenter(2);
xw = tracker_info.fixationSymbolSize(1);
yh = tracker_info.fixationSymbolSize(2);
fixation_bbox = [xc - xw/2, yc - yh/2, xc + xw/2, yc + yh/2];
if strcmpi(tracker_info.fixationSymbol, 'r')
    Screen('FillRect', wPtr, color, fixation_bbox);
elseif strcmpi(tracker_info.fixationSymbol, 'c')
    Screen('FillOval', wPtr, color, fixation_bbox);
else
    Screen('DrawLine', wPtr, color, fixation_bbox(1), yc, fixation_bbox(3), yc);
    Screen('DrawLine', wPtr, color, xc, fixation_bbox(2), xc, fixation_bbox(4));
end