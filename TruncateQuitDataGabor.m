function [GaborData, image_collection] = TruncateQuitDataGabor(GaborData, image_collection)
%TRUNCATEQUITDATAGABOR given that the subject quit the experiment before
%completion, return a truncated version of GaborData and the associated
%image_collection
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
end