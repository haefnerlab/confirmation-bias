function SaveQuitDataGabor(subjectID, phase, directory)
%% Very important to remember that you should save till BEFORE current trial and not upto current trial!!That is consider_trials<current_trial

if phase == 0
    expt_type = 'Contrast';
elseif phase == 1
    expt_type = 'Ratio';
else
    error('Unknown phase: %d', phase);
end

filename = fullfile(directory, 'RawData', sprintf('%s-GaborData%s.mat', subjectID, expt_type));
filenameQuit = fullfile(directory, 'RawData', sprintf('%s-GaborData%sQuit.mat', subjectID, expt_type));
if ~exist(filenameQuit, 'file')
    error('ERROR! Missing File: %s\nMaybe the Preliminary phase is saved under a different name?', filenameQuit);
else
    contents = load(filenameQuit);
    GaborData = contents.GaborData;
    image_collection = contents.image_collection;
end

% Truncate per-trial fields
point = GaborData.current_trial-2;
GaborData.current_trial = point;

GaborData.contrast(point+1:end) = [];
GaborData.ratio(point+1:end) = [];
GaborData.pixel_noise(point+1:end) = [];
GaborData.step_size (point+1:end) = [];

GaborData.streak (point+1:end) = [];
GaborData.reversal_counter (point+1:end) = [];
GaborData.correct_answer(point+1:end) = [];
GaborData.ideal_answer(point+1:end) = [];
GaborData.reaction_time(point+1:end) = [];
GaborData.choice(point+1:end) = [];
GaborData.accuracy(point+1:end) = [];
GaborData.order_of_orientations(point+1:end,:) = [];
GaborData.log_frame_odds(point+1:end,:) = [];
GaborData.log_decision_odds(point+1:end,:) = [];
GaborData.average_orientations(:,point+1) = [];

GaborData.eye_tracker_points(point+1:end) = [];

image_collection(point+1:end,:,:,:)=[];

save(filename, 'GaborData','image_collection');
end