function [GaborData, image_collection] = ConcatGaborData(Data1, images1, Data2, images2)
%CONCATPRELIMGABOR Concatenate experiment results.

GaborData = Data1;

GaborData.current_trial = Data1.current_trial + Data2.current_trial;

GaborData.contrast = [Data1.contrast, Data2.contrast];
GaborData.ratio = [Data1.ratio, Data2.ratio];
GaborData.pixel_noise = [Data1.pixel_noise, Data2.pixel_noise];
GaborData.step_size = [Data1.step_size, Data2.step_size];

GaborData.streak = [Data1.streak, Data2.streak];
GaborData.reversal_counter = [Data1.reversal_counter, Data2.reversal_counter ];
GaborData.correct_answer = [Data1.correct_answer, Data2.correct_answer];
GaborData.ideal_answer = [Data1.ideal_answer, Data2.ideal_answer];
GaborData.reaction_time = [Data1.reaction_time, Data2.reaction_time];
GaborData.choice = [Data1.choice, Data2.choice];
GaborData.accuracy = [Data1.accuracy, Data2.accuracy];
GaborData.order_of_orientations = [Data1.order_of_orientations; Data2.order_of_orientations];
GaborData.log_frame_odds = [Data1.log_frame_odds; Data2.log_frame_odds];
GaborData.log_decision_odds = [Data1.log_decision_odds; Data2.log_decision_odds];
GaborData.average_orientations = [Data1.average_orientations, Data2.average_orientations];

GaborData.eye_tracker_points = horzcat(Data1.eye_tracker_points, Data2.eye_tracker_points);

image_collection = cat(1, images1, images2);
end