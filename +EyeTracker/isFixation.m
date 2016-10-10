function fixation = isFixation(tracker_info, xy_pixel)
if all(tracker_info.fixationCorrection == [0 0])
    warning('Are you sure there is no correction to the fixation? (Did you remember to do [.., tracker_info, ..] = getFixation(tracker_info)?');
end
n_points = size(xy_pixel, 1);
fixation = false(n_points, 1);
target_bbox = ptbCenteredRect(tracker_info.fixationCenter - tracker_info.fixationCorrection, tracker_info.fixationRect);
for i=1:n_points
    xp = xy_pixel(i, 1);
    yp = xy_pixel(i, 2);
    fixation(i) = ...
        xp >= target_bbox(1) && yp >= target_bbox(2) && ...
        xp <= target_bbox(3) && yp <= target_bbox(4);
end
end