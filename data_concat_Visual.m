function [] = data_concat_Visual(subjectID1, subjectID2, subjectID_output, directory)

% This function will take two data files and combine them together, but it can only do two at a time
% subjectID1 and subjectID2 are the two files being combined
% subjectID_output is the desired name of the output file and can be combined with data from another
% experiment session
% directory allows this code to be able to create and save files of the subject data on any computer

% Example command line input: data_concat_Visual('MatthewFirstSession', 'MatthewSecondSession', MatthewCombinedSession', '/Users/bcs206/Documents/Summer/RawData/')

fileName1 = sprintf('%s%s-VisualTest.mat', directory, subjectID1);   %Only concat test files, never the prelims
fileName2 = sprintf('%s%s-VisualTest.mat', directory, subjectID2);
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
Test_Data.average_flashes = [Data1.Test_Data.average_flashes, Data2.Test_Data.average_flashes];
Test_Data.contrast = [Data1.Test_Data.contrast, Data2.Test_Data.contrast];
Test_Data.correct_answer = [Data1.Test_Data.correct_answer, Data2.Test_Data.correct_answer];
Test_Data.reaction_time = [Data1.Test_Data.reaction_time, Data2.Test_Data.reaction_time];
Test_Data.choice = [Data1.Test_Data.choice, Data2.Test_Data.choice];
Test_Data.accuracy = [Data1.Test_Data.accuracy, Data2.Test_Data.accuracy];

Test_Data.order_of_flashes = [Data1.Test_Data.order_of_flashes; Data2.Test_Data.order_of_flashes];  % squeeze?
			% Note this is a row-wise concatenation instead of column-wise

Test_Data.screen_frame = Data1.Test_Data.screen_frame;
Test_Data.image_length_x = Data1.Test_Data.image_length_x;
Test_Data.image_length_y = Data1.Test_Data.image_length_y;
Test_Data.screen_resolution = Data1.Test_Data.screen_resolution;

Test_Data.flash_rate = [Data1.Test_Data.flash_rate, Data2.Test_Data.flash_rate];


fileNameEnd = sprintf('%s%s-VisualTest.mat', directory, subjectID_output);
save(fileNameEnd, 'Test_Data');
end