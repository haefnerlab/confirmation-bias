function drawFixationSymbol(tracker_info, wPtr)
xc = tracker_info.fixationCenter(1);
yc = tracker_info.fixationCenter(2);
xw = tracker_info.fixationSymbolSize(1);
yh = tracker_info.fixationSymbolSize(2);
fixation_bbox = [xc - xw/2, yc - yh/2, xc + xw/2, yc + yh/2];
if strcmpi(tracker_info.fixationSymbol, 'r')
    Screen('FillRect', wPtr, tracker_info.fixationSymbolColors(1, :), fixation_bbox);
elseif strcmpi(tracker_info.fixationSymbol, 'c')
    Screen('FillOval', wPtr, tracker_info.fixationSymbolColors(1, :), fixation_bbox);
elseif strcmpi(tracker_info.fixationSymbol, 'b')
    outer_size = tracker_info.fixationSymbolSize(1);
    inner_size = round(0.3 * outer_size);
    Screen('DrawDots', wPtr, tracker_info.fixationCenter, outer_size, tracker_info.fixationSymbolColors(2, :), [0 0], 1);
    Screen('DrawDots', wPtr, tracker_info.fixationCenter, inner_size, tracker_info.fixationSymbolColors(1, :), [0 0], 1);
    Screen('DrawLine', wPtr, tracker_info.fixationSymbolColors(1, :), xc-outer_size/2, yc, xc+outer_size/2, yc, 1);
    Screen('DrawLine', wPtr, tracker_info.fixationSymbolColors(1, :), xc, yc-outer_size/2, xc, yc+outer_size/2, 1);
else
    Screen('DrawLine', wPtr, tracker_info.fixationSymbolColors(1, :), fixation_bbox(1), yc, fixation_bbox(3), yc);
    Screen('DrawLine', wPtr, tracker_info.fixationSymbolColors(1, :), xc, fixation_bbox(2), xc, fixation_bbox(4));
end