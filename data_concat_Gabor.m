function [] = data_concat_Gabor(subjectID1, subjectID2, subjectID_output, phase, directory)

% This function will take two data files and combine them together, but it can only do two at a time
% subjectID1 and subjectID2 are the two files being combined
% subjectID_output is the desired name of the output file and can be combined with data from another
% experiment session
% phase=0 means contrast concatenation and phase=1 means ratio concatenation
% directory allows this code to be able to create and save files of the subject data on any computer

% Example command line input: data_concat_Gabor('MatthewFirstSession', 'MatthewSecondSession', MatthewCombinedSession', '/Users/bcs206/Documents/Summer/RawData')

%fileName1 = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']);   %Only concat test files, never the prelims
%fileName2 = fullfile(directory, 'RawData', [subjectID '-GaborTest.mat']);

if phase==0 || phase==2
    filename1 = fullfile(directory, 'RawData', [subjectID1 '-GaborTestContrast.mat']);
    filename2 = fullfile(directory, 'RawData', [subjectID2 '-GaborTestContrast.mat']);
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Test_Data
    end
    
    Test_Data.current_trial = Data1.Test_Data.current_trial + Data2.Test_Data.current_trial;
    Test_Data.number_of_images = Data1.Test_Data.number_of_images;
    Test_Data.contrast = [Data1.Test_Data.contrast, Data2.Test_Data.contrast];
    Test_Data.correct_answer = [Data1.Test_Data.correct_answer, Data2.Test_Data.correct_answer];
    Test_Data.reaction_time = [Data1.Test_Data.reaction_time, Data2.Test_Data.reaction_time];
    Test_Data.choice = [Data1.Test_Data.choice, Data2.Test_Data.choice];
    Test_Data.accuracy = [Data1.Test_Data.accuracy, Data2.Test_Data.accuracy];
    
    Test_Data.order_of_orientations = [Data1.Test_Data.order_of_orientations; Data2.Test_Data.order_of_orientations];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    Test_Data.screen_frame = Data1.Test_Data.screen_frame;
    Test_Data.image_length_x = Data1.Test_Data.image_length_x;
    Test_Data.image_length_y = Data1.Test_Data.image_length_y;
    Test_Data.screen_resolution = Data1.Test_Data.screen_resolution;
    
    Test_Data.log_odds = [Data1.Test_Data.log_odds, Data2.Test_Data.log_odds];
    Test_Data.ratios = [Data1.Test_Data.ratios, Data2.Test_Data.ratios];
    
    Test_Data.image_template1 = [Data1.Test_Data.image_template1; Data2.Test_Data.image_template1];  % squeeze?
    Test_Data.image_template2 = [Data1.Test_Data.template2; Data2.Test_Data.template2];  % squeeze?
    Test_Data.image_template_difference = [Data1.Test_Data.image_template_difference; Data2.Test_Data.image_template_difference];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    
    [trials1, frames1, length_x1, length_y1] = size(Data1.image_collection);
    [trials2, frames2, length_x2, length_y2] = size(Data2.image_collection);
    image_collection = zeros(trials1+trials2, frames1, length_x1, length_y1);
    image_collection(1:trials1, :, :, :) = Data1.image_collection;
    image_collection(trials1+1:end, :, :, :) = Data2.image_collection;
    
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID '-GaborTestContrast.mat']);
    save(fileNameEnd, 'image_collection', 'Test_Data');
end

if phase==1 || phase==2
    filename1 = fullfile(directory, 'RawData', [subjectID1 '-GaborTestRatio.mat']);
    filename2 = fullfile(directory, 'RawData', [subjectID2 '-GaborTestRatio.mat']);
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Test_Data
    end
    
    Test_Data.current_trial = Data1.Test_Data.current_trial + Data2.Test_Data.current_trial;
    Test_Data.number_of_images = Data1.Test_Data.number_of_images;
    Test_Data.contrast = [Data1.Test_Data.contrast, Data2.Test_Data.contrast];
    Test_Data.correct_answer = [Data1.Test_Data.correct_answer, Data2.Test_Data.correct_answer];
    Test_Data.reaction_time = [Data1.Test_Data.reaction_time, Data2.Test_Data.reaction_time];
    Test_Data.choice = [Data1.Test_Data.choice, Data2.Test_Data.choice];
    Test_Data.accuracy = [Data1.Test_Data.accuracy, Data2.Test_Data.accuracy];
    
    Test_Data.order_of_orientations = [Data1.Test_Data.order_of_orientations; Data2.Test_Data.order_of_orientations];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    Test_Data.screen_frame = Data1.Test_Data.screen_frame;
    Test_Data.image_length_x = Data1.Test_Data.image_length_x;
    Test_Data.image_length_y = Data1.Test_Data.image_length_y;
    Test_Data.screen_resolution = Data1.Test_Data.screen_resolution;
    
    Test_Data.log_odds = [Data1.Test_Data.log_odds, Data2.Test_Data.log_odds];
    Test_Data.ratios = [Data1.Test_Data.ratios, Data2.Test_Data.ratios];
    
    Test_Data.image_template1 = [Data1.Test_Data.image_template1; Data2.Test_Data.image_template1];  % squeeze?
    Test_Data.image_template2 = [Data1.Test_Data.template2; Data2.Test_Data.template2];  % squeeze?
    Test_Data.image_template_difference = [Data1.Test_Data.image_template_difference; Data2.Test_Data.image_template_difference];  % squeeze?
    % Note this is a row-wise concatenation instead of column-wise
    
    
    [trials1, frames1, length_x1, length_y1] = size(Data1.image_collection);
    [trials2, frames2, length_x2, length_y2] = size(Data2.image_collection);
    image_collection = zeros(trials1+trials2, frames1, length_x1, length_y1);
    image_collection(1:trials1, :, :, :) = Data1.image_collection;
    image_collection(trials1+1:end, :, :, :) = Data2.image_collection;
    
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID '-GaborTestRatio.mat']);
    save(fileNameEnd, 'image_collection', 'Test_Data');
end
end