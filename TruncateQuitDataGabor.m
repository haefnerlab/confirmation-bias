function GaborData = TruncateQuitDataGabor(GaborData)
%TRUNCATEQUITDATAGABOR given that the subject quit the experiment before
%completion, return a truncated version of GaborData
point = GaborData.current_trial-2;
GaborData.current_trial = point;

GaborData.contrast(point+1:end) = [];
GaborData.ratio(point+1:end) = [];
GaborData.noise(point+1:end) = [];
GaborData.step_size (point+1:end) = [];

GaborData.iid(point+1:end) = [];
GaborData.seed(point+1:end) = [];
if isfield(GaborData, 'checksum'), GaborData.checksum(point+1:end) = []; end
GaborData.streak (point+1:end) = [];
GaborData.reversal_counter (point+1:end) = [];
GaborData.correct_answer(point+1:end) = [];
GaborData.ideal_answer(point+1:end) = [];
GaborData.reaction_time(point+1:end) = [];
GaborData.choice(point+1:end) = [];
GaborData.accuracy(point+1:end) = [];
GaborData.frame_categories(point+1:end,:) = [];
GaborData.ideal_frame_signals(point+1:end,:) = [];

GaborData.eye_tracker_points(point+1:end) = [];
end