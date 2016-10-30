function [] = data_concat_Auditory(subjectID1, subjectID2, subjectID_output, phase, directory)

% This function will take two data files and combine them together, but it can only do two at a time

% subjectID1 and subjectID2 are the two files being combined

% subjectID_output is the desired name of the output file and can be combined with data from another
% experiment session

% phase determines if we are concatenating data for the volume or ratio
% test phase only, or both (0 is volume only, 1 is ratio only,2 is both )

% directory allows this code to be able to create and save files of the subject data on any computer

% Example command line input: data_concat_Audtiory('MatthewFirstSession', 'MatthewSecondSession', MatthewCombinedSession', '/Users/bcs206/Documents/Summer/RawData')

% Concatenate the volume test phases

if phase == 0 || phase == 2
    fileName1 = fullfile(directory, 'RawData', [subjectID1 '-AuditoryTestVolume.mat']);
    fileName2 = fullfile(directory, 'RawData', [subjectID2 '-AuditoryTestVolume.mat']);
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Test_Data
    end
    
    Test_Data.move_on = [Data1.Test_Data.move_on, Data2.Test_Data.move_on];
    Test_Data.step_size = [Data1.Test_Data.step_size, Data2.Test_Data.step_size];
    Test_Data.reversal_counter = [Data1.Test_Data.reversal_counter, Data2.Test_Data.reversal_counter ];
    Test_Data.volume = [Data1.Test_Data.volume, Data2.Test_Data.volume];
    Test_Data.correct_answer = [Data1.Test_Data.correct_answer, Data2.Test_Data.correct_answer];
    Test_Data.staircase_answer = [Data1.Test_Data.staircase_answer, Data2.Test_Data.staircase_answer];
    Test_Data.reaction_time = [Data1.Test_Data.reaction_time, Data2.Test_Data.reaction_time];
    Test_Data.choice = [Data1.Test_Data.choice, Data2.Test_Data.choice];
    Test_Data.accuracy = [Data1.Test_Data.accuracy, Data2.Test_Data.accuracy];
    Test_Data.ratio = [Data1.Test_Data.ratio, Data2.Test_Data.ratio];
    Test_Data.test_phase = [Data1.Test_Data.test_phase, Data2.Test_Data.test_phase];
    
    Test_Data.current_trial = Data1.Test_Data.current_trial + Data2.Test_Data.current_trial;
    Test_Data.sampling_rate = Data1.Test_Data.sampling_rate;
    Test_Data.stimulus_duration = Data1.Test_Data.stimulus_duration;
    Test_Data.bins = Data1.Test_Data.bins;
    Test_Data.number_of_frames = Data1.Test_Data.number_of_frames;
    
    Test_Data.click_rate = [Data1.Test_Data.click_rate, Data2.Test_Data.click_rate];
    Test_Data.average_clicks = [Data1.Test_Data.average_clicks, Data2.Test_Data.average_clicks];
    
    Test_Data.order_of_clicks = [Data1.Test_Data.order_of_clicks; Data2.Test_Data.order_of_clicks];
    % Note this is a row-wise concatenation instead of column-wise
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID '-AuditoryTestVolume.mat']);
    save(fileNameEnd, 'Test_Data');
end    
if phase == 1 || phase == 2
    % Now concatenate the ratio test phases
    
    fileName1 = fullfile(directory, 'RawData', [subjectID1 '-AuditoryTestRatio.mat']);
    fileName2 = fullfile(directory, 'RawData', [subjectID2 '-AuditoryTestRatio.mat']);
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Test_Data
    end
    
    Test_Data.move_on = [Data1.Test_Data.move_on, Data2.Test_Data.move_on];
    Test_Data.step_size = [Data1.Test_Data.step_size, Data2.Test_Data.step_size];
    Test_Data.reversal_counter = [Data1.Test_Data.reversal_counter, Data2.Test_Data.reversal_counter ];
    Test_Data.volume = [Data1.Test_Data.volume, Data2.Test_Data.volume];
    Test_Data.correct_answer = [Data1.Test_Data.correct_answer, Data2.Test_Data.correct_answer];
    Test_Data.staircase_answer = [Data1.Test_Data.staircase_answer, Data2.Test_Data.staircase_answer];
    Test_Data.reaction_time = [Data1.Test_Data.reaction_time, Data2.Test_Data.reaction_time];
    Test_Data.choice = [Data1.Test_Data.choice, Data2.Test_Data.choice];
    Test_Data.accuracy = [Data1.Test_Data.accuracy, Data2.Test_Data.accuracy];
    Test_Data.ratio = [Data1.Test_Data.ratio, Data2.Test_Data.ratio];
    Test_Data.test_phase = [Data1.Test_Data.test_phase, Data2.Test_Data.test_phase];
    
    Test_Data.current_trial = Data1.Test_Data.current_trial + Data2.Test_Data.current_trial;
    Test_Data.sampling_rate = Data1.Test_Data.sampling_rate;
    Test_Data.stimulus_duration = Data1.Test_Data.stimulus_duration;
    Test_Data.bins = Data1.Test_Data.bins;
    Test_Data.number_of_frames = Data1.Test_Data.number_of_frames;
    
    Test_Data.click_rate = [Data1.Test_Data.click_rate, Data2.Test_Data.click_rate];
    Test_Data.average_clicks = [Data1.Test_Data.average_clicks, Data2.Test_Data.average_clicks];
    
    Test_Data.order_of_clicks = [Data1.Test_Data.order_of_clicks; Data2.Test_Data.order_of_clicks];
    % Note this is a row-wise concatenation instead of column-wise
    
    
    fileNameEnd = fullfile(directory, 'RawData', [subjectID '-AuditoryTestRatio.mat']);
    save(fileNameEnd, 'Test_Data');
    
end
end