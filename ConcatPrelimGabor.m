function [] = ConcatPrelimGabor(subjectID1, subjectID2, subjectID_output, phase, directory)

% This function will take two data files and combine them together, but it can only do two at a time
% subjectID1 and subjectID2 are the two files being combined
% subjectID_output is the desired name of the output file and can be combined with data from another
% experiment session
% phase=0 means contrast concatenation and phase=1 means ratio concatenation
% directory allows this code to be able to create and save files of the subject data on any computer

% Example command line input: data_concat_Gabor('MatthewFirstSession', 'MatthewSecondSession', MatthewCombinedSession', '/Users/bcs206/Documents/Summer/RawData')

%fileName1 = fullfile(directory, 'RawData', [subjectID '-GaborPreliminary.mat']);   %Only concat Preliminary files, never the prelims
%fileName2 = fullfile(directory, 'RawData', [subjectID '-GaborPreliminary.mat']);

if phase==0 || phase==2
    fileName1 = fullfile(directory, 'RawData', [subjectID1 '-GaborDataContrast.mat']);
    fileName2 = fullfile(directory, 'RawData', [subjectID2 '-GaborDataContrast.mat']);
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Preliminary_Data
    end
    
    Preliminary_Data.current_trial = Data1.Preliminary_Data.current_trial + Data2.Preliminary_Data.current_trial;
    Preliminary_Data.number_of_images = Data1.Preliminary_Data.number_of_images;
    Preliminary_Data.contrast = [Data1.Preliminary_Data.contrast, Data2.Preliminary_Data.contrast];
    Preliminary_Data.correct_answer = [Data1.Preliminary_Data.correct_answer, Data2.Preliminary_Data.correct_answer];
    Preliminary_Data.reaction_time = [Data1.Preliminary_Data.reaction_time, Data2.Preliminary_Data.reaction_time];
    Preliminary_Data.choice = [Data1.Preliminary_Data.choice, Data2.Preliminary_Data.choice];
    Preliminary_Data.accuracy = [Data1.Preliminary_Data.accuracy, Data2.Preliminary_Data.accuracy];
    
    Preliminary_Data.order_of_orientations = [Data1.Preliminary_Data.order_of_orientations; Data2.Preliminary_Data.order_of_orientations];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    Preliminary_Data.screen_frame = Data1.Preliminary_Data.screen_frame;
    Preliminary_Data.image_length_x = Data1.Preliminary_Data.image_length_x;
    Preliminary_Data.image_length_y = Data1.Preliminary_Data.image_length_y;
    Preliminary_Data.screen_resolution = Data1.Preliminary_Data.screen_resolution;
    
    Preliminary_Data.log_odds = [Data1.Preliminary_Data.log_odds, Data2.Preliminary_Data.log_odds];
    Preliminary_Data.ratio = [Data1.Preliminary_Data.ratio, Data2.Preliminary_Data.ratio];
    
    Preliminary_Data.image_template1 = [Data1.Preliminary_Data.image_template1; Data2.Preliminary_Data.image_template1];  % squeeze?
    Preliminary_Data.image_template2 = [Data1.Preliminary_Data.image_template2; Data2.Preliminary_Data.image_template2];  % squeeze?
    Preliminary_Data.image_template_difference = [Data1.Preliminary_Data.image_template_difference; Data2.Preliminary_Data.image_template_difference];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    
    Preliminary_Data.move_on = [Data1.Preliminary_Data.move_on, Data2.Preliminary_Data.move_on];
    Preliminary_Data.step_size = [Data1.Preliminary_Data.step_size, Data2.Preliminary_Data.step_size];
    Preliminary_Data.reversal_counter = [Data1.Preliminary_Data.reversal_counter, Data2.Preliminary_Data.reversal_counter ];
    Preliminary_Data.staircase_answer = [Data1.Preliminary_Data.staircase_answer, Data2.Preliminary_Data.staircase_answer];
    Preliminary_Data.average_orientations = [Data1.Preliminary_Data.average_orientations, Data2.Preliminary_Data.average_orientations];
    %%Preliminary_Data.eye_tracker_points = {};
    Preliminary_Data.left_template = Data1.Preliminary_Data.left_template;
    Preliminary_Data.right_template = Data1.Preliminary_Data.right_template;
  
    
    [trials1, frames1, length_x1, length_y1] = size(Data1.image_collection);
    [trials2, frames2, length_x2, length_y2] = size(Data2.image_collection);
    image_collection = zeros(trials1+trials2, frames1, length_x1, length_y1);
    image_collection(1:trials1, :, :, :) = Data1.image_collection;
    image_collection(trials1+1:end, :, :, :) = Data2.image_collection;
    
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID_output '-GaborDataContrast.mat']);
    save(fileNameEnd, 'image_collection', 'Preliminary_Data');
end

if phase==1 || phase==2
    fileName1 = fullfile(directory, 'RawData', [subjectID1 '-GaborDataRatio.mat']);
    fileName2 = fullfile(directory, 'RawData', [subjectID2 '-GaborDataRatio.mat']);
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Preliminary_Data
    end
    
    Preliminary_Data.current_trial = Data1.Preliminary_Data.current_trial + Data2.Preliminary_Data.current_trial;
    Preliminary_Data.number_of_images = Data1.Preliminary_Data.number_of_images;
    Preliminary_Data.contrast = [Data1.Preliminary_Data.contrast, Data2.Preliminary_Data.contrast];
    Preliminary_Data.correct_answer = [Data1.Preliminary_Data.correct_answer, Data2.Preliminary_Data.correct_answer];
    Preliminary_Data.reaction_time = [Data1.Preliminary_Data.reaction_time, Data2.Preliminary_Data.reaction_time];
    Preliminary_Data.choice = [Data1.Preliminary_Data.choice, Data2.Preliminary_Data.choice];
    Preliminary_Data.accuracy = [Data1.Preliminary_Data.accuracy, Data2.Preliminary_Data.accuracy];
    
    Preliminary_Data.order_of_orientations = [Data1.Preliminary_Data.order_of_orientations; Data2.Preliminary_Data.order_of_orientations];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    Preliminary_Data.screen_frame = Data1.Preliminary_Data.screen_frame;
    Preliminary_Data.image_length_x = Data1.Preliminary_Data.image_length_x;
    Preliminary_Data.image_length_y = Data1.Preliminary_Data.image_length_y;
    Preliminary_Data.screen_resolution = Data1.Preliminary_Data.screen_resolution;
    
    Preliminary_Data.log_odds = [Data1.Preliminary_Data.log_odds, Data2.Preliminary_Data.log_odds];
    Preliminary_Data.ratio = [Data1.Preliminary_Data.ratio, Data2.Preliminary_Data.ratio];
    
    Preliminary_Data.image_template1 = [Data1.Preliminary_Data.image_template1; Data2.Preliminary_Data.image_template1];  % squeeze?
    Preliminary_Data.image_template2 = [Data1.Preliminary_Data.template2; Data2.Preliminary_Data.template2];  % squeeze?
    Preliminary_Data.image_template_difference = [Data1.Preliminary_Data.image_template_difference; Data2.Preliminary_Data.image_template_difference];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    Preliminary_Data.move_on = [Data1.Preliminary_Data.move_on, Data2.Preliminary_Data.move_on];
    Preliminary_Data.step_size = [Data1.Preliminary_Data.step_size, Data2.Preliminary_Data.step_size];
    Preliminary_Data.reversal_counter = [Data1.Preliminary_Data.reversal_counter, Data2.Preliminary_Data.reversal_counter ];
    Preliminary_Data.staircase_answer = [Data1.Preliminary_Data.staircase_answer, Data2.Preliminary_Data.staircase_answer];
    Preliminary_Data.average_orientations = [Data1.Preliminary_Data.average_orientations, Data2.Preliminary_Data.average_orientations];
    %%Preliminary_Data.eye_tracker_points = {};
    Preliminary_Data.left_template = Data1.Preliminary_Data.left_template;
    Preliminary_Data.right_template = Data1.Preliminary_Data.right_template;
  
    
    
    [trials1, frames1, length_x1, length_y1] = size(Data1.image_collection);
    [trials2, frames2, length_x2, length_y2] = size(Data2.image_collection);
    image_collection = zeros(trials1+trials2, frames1, length_x1, length_y1);
    image_collection(1:trials1, :, :, :) = Data1.image_collection;
    image_collection(trials1+1:end, :, :, :) = Data2.image_collection;
    
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID_output '-GaborDataRatio.mat']);
    save(fileNameEnd, 'image_collection', 'Preliminary_Data');
end
end