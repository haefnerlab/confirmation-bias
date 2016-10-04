function is_valid = isValidTrial(tracker_info, trial_data)
is_valid = EyeTracker.didHoldFixation(tracker_info, trial_data) && ...
    ~EyeTracker.didBlink(tracker_info, trial_data);
end