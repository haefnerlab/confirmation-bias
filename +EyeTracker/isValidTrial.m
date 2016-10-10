function is_valid = isValidTrial(tracker_info, trial_xy_data)
is_valid = all(EyeTracker.isFixation(tracker_info, trial_xy_data));
end