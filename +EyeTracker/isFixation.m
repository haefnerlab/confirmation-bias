function fixation = isFixation(tracker_info, xy_pixel)
if all(tracker_info.fixationCorrection == [0 0])
    warning('Are you sure there is no correction to the fixation? (Did you remember to do [.., tracker_info, ..] = getFixation(tracker_info)?');
end
n_points = size(xy_pixel, 1);
xy_adjusted = xy_pixel + tracker_info.fixationCorrection;
diffs = xy_adjusted - repmat(tracker_info.fixationCenter, n_points, 1);
distances = sqrt(sum(diffs.^2, 2));
fixation = distances < tracker_info.fixationRadius;
end