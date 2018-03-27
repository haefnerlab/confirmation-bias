function GaborData = ConcatGaborData(Data1, Data2)
%CONCATPRELIMGABOR Concatenate experiment results.

GaborData = Data1;

GaborData.current_trial = Data1.current_trial + Data2.current_trial;

GaborData.contrast = [Data1.contrast, Data2.contrast];
GaborData.ratio = [Data1.ratio, Data2.ratio];
GaborData.noise = [Data1.noise, Data2.noise];
GaborData.step_size = [Data1.step_size, Data2.step_size];

GaborData.iid = [Data1.iid, Data2.iid];
GaborData.seed = [Data1.seed, Data2.seed];
if isfield(GaborData, 'checksum'), GaborData.checksum = [Data1.checksum, Data2.checksum]; end
GaborData.streak = [Data1.streak, Data2.streak];
GaborData.reversal_counter = [Data1.reversal_counter, Data2.reversal_counter ];
GaborData.correct_answer = [Data1.correct_answer, Data2.correct_answer];
GaborData.ideal_answer = [Data1.ideal_answer, Data2.ideal_answer];
GaborData.reaction_time = [Data1.reaction_time, Data2.reaction_time];
GaborData.choice = [Data1.choice, Data2.choice];
GaborData.accuracy = [Data1.accuracy, Data2.accuracy];
GaborData.frame_categories = [Data1.frame_categories; Data2.frame_categories];
GaborData.ideal_frame_signals = [Data1.ideal_frame_signals; Data2.ideal_frame_signals];

GaborData.eye_tracker_points = horzcat(Data1.eye_tracker_points, Data2.eye_tracker_points);
end
