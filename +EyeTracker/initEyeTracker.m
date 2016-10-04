function tracker_info = initEyeTracker()
% Connect to eye tracker.
% Return struct with tracker's configuration options.
tracker_info = struct(...
    'coordinates', 'pixel');
end