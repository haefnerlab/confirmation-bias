function ConcatPrelimGabor(subjectID1, subjectID2, subjectID_output, phase, directory)
% This function will take two data files and combine them together, but it can only do two at a time
% subjectID1 and subjectID2 are the two files being combined
% subjectID_output is the desired name of the output file and can be combined with data from another
% experiment session
% phase=0 means contrast concatenation and phase=1 means ratio concatenation
% directory allows this code to be able to create and save files of the subject data on any computer

% Example command line input: data_concat_Gabor('MatthewFirstSession', 'MatthewSecondSession', MatthewCombinedSession', '/Users/bcs206/Documents/Summer/RawData')

%fileName1 = fullfile(directory, 'RawData', [subjectID '-GaborPreliminary.mat']);   %Only concat Preliminary files, never the prelims
%fileName2 = fullfile(directory, 'RawData', [subjectID '-GaborPreliminary.mat']);

if phase == 0
    expt_type = 'Contrast';
elseif phase == 1
    expt_type = 'Ratio';
else
    error('Unknown phase: %d', phase);
end

fileName1 = fullfile(directory, 'RawData', sprintf('%s-GaborData%s.mat', subjectID1, expt_type));
fileName2 = fullfile(directory, 'RawData', sprintf('%s-GaborData%s.mat', subjectID2, expt_type));

if ~exist(fileName1, 'file')
    error('Missing File: %s', fileName1);
elseif ~exist(fileName2, 'file')
    error('Missing File: %s', fileName2);
else
    Data1 = load(fileName1);
    Data2 = load(fileName2);
    % This will save the data to two different structs, each containing
    % a different GaborData struct
end

GaborData = Data1.GaborData;

GaborData.current_trial = Data1.GaborData.current_trial + Data2.GaborData.current_trial;

GaborData.contrast = [Data1.GaborData.contrast, Data2.GaborData.contrast];
GaborData.ratio = [Data1.GaborData.ratio, Data2.GaborData.ratio];
GaborData.pixel_noise = [Data1.GaborData.pixel_noise, Data2.GaborData.pixel_noise];
GaborData.step_size = [Data1.GaborData.step_size, Data2.GaborData.step_size];

GaborData.streak = [Data1.GaborData.streak, Data2.GaborData.streak];
GaborData.reversal_counter = [Data1.GaborData.reversal_counter, Data2.GaborData.reversal_counter ];
GaborData.correct_answer = [Data1.GaborData.correct_answer, Data2.GaborData.correct_answer];
GaborData.ideal_answer = [Data1.GaborData.ideal_answer, Data2.GaborData.ideal_answer];
GaborData.reaction_time = [Data1.GaborData.reaction_time, Data2.GaborData.reaction_time];
GaborData.choice = [Data1.GaborData.choice, Data2.GaborData.choice];
GaborData.accuracy = [Data1.GaborData.accuracy, Data2.GaborData.accuracy];
GaborData.order_of_orientations = [Data1.GaborData.order_of_orientations; Data2.GaborData.order_of_orientations];
GaborData.log_frame_odds = [Data1.GaborData.log_frame_odds; Data2.GaborData.log_frame_odds];
GaborData.log_decision_odds = [Data1.GaborData.log_decision_odds; Data2.GaborData.log_decision_odds];
GaborData.average_orientations = [Data1.GaborData.average_orientations, Data2.GaborData.average_orientations];

GaborData.eye_tracker_points = horzcat(Data1.GaborData.eye_tracker_points, Data2.GaborData.eye_tracker_points);

image_collection = cat(1, Data1.image_collection, Data2.image_collection);

fileNameEnd = fullfile(directory, 'RawData', sprintf('%s-GaborData%s.mat', subjectID_output, expt_type));
save(fileNameEnd, 'image_collection', 'GaborData');
end