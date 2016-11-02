function [] = ConcatPrelimAuditory(subjectID1, subjectID2, subjectID_output,phase, directory)
if phase==0
    fileName1 = [directory 'RawData/' subjectID1 '-AuditoryDataVolume.mat'];
    fileName2 = [directory 'RawData/' subjectID2 '-AuditoryDataVolume.mat'];
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Preliminary_Data
    end
    Preliminary_Data.move_on = [Data1.Preliminary_Data.move_on, Data2.Preliminary_Data.move_on];
    Preliminary_Data.step_size = [Data1.Preliminary_Data.step_size, Data2.Preliminary_Data.step_size];
    Preliminary_Data.reversal_counter = [Data1.Preliminary_Data.reversal_counter, Data2.Preliminary_Data.reversal_counter ];
    Preliminary_Data.volume = [Data1.Preliminary_Data.volume, Data2.Preliminary_Data.volume];
    Preliminary_Data.correct_answer = [Data1.Preliminary_Data.correct_answer, Data2.Preliminary_Data.correct_answer];
    Preliminary_Data.staircase_answer = [Data1.Preliminary_Data.staircase_answer, Data2.Preliminary_Data.staircase_answer];
    Preliminary_Data.reaction_time = [Data1.Preliminary_Data.reaction_time, Data2.Preliminary_Data.reaction_time];
    Preliminary_Data.choice = [Data1.Preliminary_Data.choice, Data2.Preliminary_Data.choice];
    Preliminary_Data.accuracy = [Data1.Preliminary_Data.accuracy, Data2.Preliminary_Data.accuracy];
    Preliminary_Data.ratio = [Data1.Preliminary_Data.ratio, Data2.Preliminary_Data.ratio];
    Preliminary_Data.test_phase = [Data1.Preliminary_Data.test_phase, Data2.Preliminary_Data.test_phase];
    
    Preliminary_Data.current_trial = Data1.Preliminary_Data.current_trial + Data2.Preliminary_Data.current_trial;
    Preliminary_Data.sampling_rate = Data1.Preliminary_Data.sampling_rate;
    Preliminary_Data.stimulus_duration = Data1.Preliminary_Data.stimulus_duration;
    Preliminary_Data.bins = Data1.Preliminary_Data.bins;
    Preliminary_Data.number_of_frames = Data1.Preliminary_Data.number_of_frames;
    
    Preliminary_Data.click_rate = [Data1.Preliminary_Data.click_rate, Data2.Preliminary_Data.click_rate];
    Preliminary_Data.average_clicks = [Data1.Preliminary_Data.average_clicks, Data2.Preliminary_Data.average_clicks];
    
    Preliminary_Data.order_of_clicks = [Data1.Preliminary_Data.order_of_clicks; Data2.Preliminary_Data.order_of_clicks];
    % Note this is a row-wise concatenation instead of column-wise
    
    
    fileNameEnd = sprintf('%s%s-AuditoryDataVolume.mat', [ directory 'RawData/'], subjectID_output);
    save(fileNameEnd, 'Preliminary_Data');
    
elseif phase == 1
    fileName1 = [directory 'RawData/' subjectID1 '-AuditoryDataRatio.mat'];
    fileName2 = [directory 'RawData/' subjectID2 '-AuditoryDataRatio.mat'];
    
    if ~exist(fileName1, 'file') || ~exist(fileName2, 'file')
        disp(strcat('ERROR! Missing File: ', fileName1));
        disp(strcat('ERROR! Missing File: ', fileName2));
        return;
    else
        Data1 = load(fileName1);
        Data2 = load(fileName2);
        % This will save the data to two different structs, each with a different Preliminary_Data
    end
    
    Preliminary_Data.move_on = [Data1.Preliminary_Data.move_on, Data2.Preliminary_Data.move_on];
    Preliminary_Data.step_size = [Data1.Preliminary_Data.step_size, Data2.Preliminary_Data.step_size];
    Preliminary_Data.reversal_counter = [Data1.Preliminary_Data.reversal_counter, Data2.Preliminary_Data.reversal_counter ];
    Preliminary_Data.volume = [Data1.Preliminary_Data.volume, Data2.Preliminary_Data.volume];
    Preliminary_Data.correct_answer = [Data1.Preliminary_Data.correct_answer, Data2.Preliminary_Data.correct_answer];
    Preliminary_Data.staircase_answer = [Data1.Preliminary_Data.staircase_answer, Data2.Preliminary_Data.staircase_answer];
    Preliminary_Data.reaction_time = [Data1.Preliminary_Data.reaction_time, Data2.Preliminary_Data.reaction_time];
    Preliminary_Data.choice = [Data1.Preliminary_Data.choice, Data2.Preliminary_Data.choice];
    Preliminary_Data.accuracy = [Data1.Preliminary_Data.accuracy, Data2.Preliminary_Data.accuracy];
    Preliminary_Data.ratio = [Data1.Preliminary_Data.ratio, Data2.Preliminary_Data.ratio];
    Preliminary_Data.test_phase = [Data1.Preliminary_Data.test_phase, Data2.Preliminary_Data.test_phase];
    
    Preliminary_Data.current_trial = Data1.Preliminary_Data.current_trial + Data2.Preliminary_Data.current_trial;
    Preliminary_Data.sampling_rate = Data1.Preliminary_Data.sampling_rate;
    Preliminary_Data.stimulus_duration = Data1.Preliminary_Data.stimulus_duration;
    Preliminary_Data.bins = Data1.Preliminary_Data.bins;
    Preliminary_Data.number_of_frames = Data1.Preliminary_Data.number_of_frames;
    
    Preliminary_Data.click_rate = [Data1.Preliminary_Data.click_rate, Data2.Preliminary_Data.click_rate];
    Preliminary_Data.average_clicks = [Data1.Preliminary_Data.average_clicks, Data2.Preliminary_Data.average_clicks];
    
    Preliminary_Data.order_of_clicks = [Data1.Preliminary_Data.order_of_clicks; Data2.Preliminary_Data.order_of_clicks];
    % Note this is a row-wise concatenation instead of column-wise
    
    
    fileNameEnd = sprintf('%s%s-AuditoryDataRatio.mat', [ directory 'RawData/'], subjectID_output);
    save(fileNameEnd, 'Preliminary_Data');
    
end


end