function [] = SaveGaborTestPhaseBelowThreshold(subjectID, phase, loops, preliminary_trials,threshold, directory)

cd(fullfile(directory, 'Code'))
if phase == 0
        filename = fullfile(directory, 'RawData', [subjectID '-GaborDataContrast.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
            return;
        else
            load(filename); % Load Preliminary_Data
        end
        
elseif phase == 1
        filename = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
            return;
        else
            load(filename); % Load Preliminary_Data
        end
end
elements=0;
if phase ==0
    for t = 1:loops*preliminary_trials
        if (Preliminary_Data.contrast(t)<=threshold)
            elements=elements+1;
        end
    end
elseif phase ==1
    for t = 1:loops*preliminary_trials
        if (Preliminary_Data.ratio(t)<=threshold)
            elements=elements+1;
        end
    end
end




Test_Data.move_on = zeros(1,elements);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
Test_Data.step_size = zeros(1,elements);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
Test_Data.reversal_counter = zeros(1,elements);   % How many trials has the subject got wrong? When to change the step size?
Test_Data.contrast = zeros(1,elements);         % How loud the sound is, or the signal level
Test_Data.number_of_images = 10;
Test_Data.correct_answer = zeros(1,elements);         % What was the right ear/answer?
Test_Data.staircase_answer = zeros(1,elements);         % What was the right ear/answer?
Test_Data.reaction_time = zeros(1,elements);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
Test_Data.choice = zeros(1,elements);                 % Which ear did the subject choose? 1 = left, 0 = right
Test_Data.accuracy = zeros(1,elements);               % Did the subject make the right selection? 1 = yes, 0 = no
Test_Data.order_of_orientations = zeros(elements,Test_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
Test_Data.log_odds = zeros(1,elements);
Test_Data.ratio = zeros(1,elements);
Test_Data.average_orientations = zeros(2,elements);

Test_Data.current_trial = 0;


Test_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
Test_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
Test_Data.image_length_x = 5;  % Size of the image along x-axis
Test_Data.image_length_y = 5;

Test_Data.image_template1 = zeros(elements,Test_Data.number_of_images);
Test_Data.image_template2 = zeros(elements,Test_Data.number_of_images);
Test_Data.image_template_difference = zeros(elements,Test_Data.number_of_images);

Test_Data.eye_tracker_points = {};

res = Test_Data.screen_resolution;
Test_Data.left_template = zeros(res * 5);
for i=1:5
    Test_Data.left_template((i-1)*res+1:i*res, (i-1)*res+1:i*res) = 1;
end
Test_Data.right_template = rot90(Test_Data.left_template);


image_collection_test = zeros(elements, Test_Data.number_of_images, ...
    Test_Data.image_length_x*Test_Data.screen_resolution, Test_Data.image_length_y*Test_Data.screen_resolution);

if phase==0
    k=1;
    for i = 1:loops*preliminary_trials
        if (Preliminary_Data.contrast(i)<=threshold)
            
            
            Test_Data.move_on(k) = Preliminary_Data.move_on(i);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
            Test_Data.step_size(k) = Preliminary_Data.step_size(i);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
            Test_Data.reversal_counter(k) = Preliminary_Data.reversal_counter(i);   % How many trials has the subject got wrong? When to change the step size?
            Test_Data.contrast(k) = Preliminary_Data.contrast(i);         % How loud the sound is, or the signal level
            Test_Data.number_of_images = 10;
            Test_Data.correct_answer(k) = Preliminary_Data.correct_answer(i);       % What was the right ear/answer?
            Test_Data.staircase_answer(k) = Preliminary_Data.staircase_answer(i);         % What was the right ear/answer?
            Test_Data.reaction_time(k) = Preliminary_Data.reaction_time(i);           % How quick is the subject to input their answer for the orientation choice? Recorded in ms
            Test_Data.choice(k) = Preliminary_Data.choice(i);                 % Which ear did the subject choose? 1 = left, 0 = right
            Test_Data.accuracy(k) = Preliminary_Data.accuracy(i);                % Did the subject make the right selection? 1 = yes, 0 = no
            Test_Data.order_of_orientations(k,:) = Preliminary_Data.order_of_orientations(i,:);         % Record the randomized order of image orientations throughout the trial
            Test_Data.log_odds(k) = Preliminary_Data.log_odds(i);
            Test_Data.ratio(k) = Preliminary_Data.ratio(i)
            Test_Data.average_orientations(:,k) = Preliminary_Data.average_orientations(:,i);
            
            Test_Data.current_trial = k;
            
            
            Test_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
            Test_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
            Test_Data.image_length_x = 5;  % Size of the image along x-axis
            Test_Data.image_length_y = 5;
            
            Test_Data.image_template1(k,:) = Preliminary_Data.image_template1(i,:);
            Test_Data.image_template2(k,:) = Preliminary_Data.image_template2(i,:);
            Test_Data.image_template_difference(k,:) = Preliminary_Data.image_template_difference(i,:);
            
            Test_Data.eye_tracker_points{k} = Preliminary_Data.eye_tracker_points{i};
            
            res = Test_Data.screen_resolution;
            Test_Data.left_template = zeros(res * 5);
            for t=1:5
                Test_Data.left_template((t-1)*res+1:t*res, (t-1)*res+1:t*res) = 1;
            end
            Test_Data.right_template = rot90(Test_Data.left_template);
            
            
            image_collection_test(k,:,:,:) = image_collection(i,:,:,:);
                
        k=k+1;
        end
    end
elseif phase==1
    k=1;
    for i = 1:loops*preliminary_trials
        if (Preliminary_Data.ratio(i)<=threshold)
            Test_Data.move_on(k) = Preliminary_Data.move_on(i);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
            Test_Data.step_size(k) = Preliminary_Data.step_size(i);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
            Test_Data.reversal_counter(k) = Preliminary_Data.reversal_counter(i);   % How many trials has the subject got wrong? When to change the step size?
            Test_Data.contrast(k) = Preliminary_Data.contrast(i);         % How loud the sound is, or the signal level
            Test_Data.number_of_images = 10;
            Test_Data.correct_answer(k) = Preliminary_Data.correct_answer(i);       % What was the right ear/answer?
            Test_Data.staircase_answer(k) = Preliminary_Data.staircase_answer(i);         % What was the right ear/answer?
            Test_Data.reaction_time(k) = Preliminary_Data.reaction_time(i);           % How quick is the subject to input their answer for the orientation choice? Recorded in ms
            Test_Data.choice(k) = Preliminary_Data.choice(i);                 % Which ear did the subject choose? 1 = left, 0 = right
            Test_Data.accuracy(k) = Preliminary_Data.accuracy(i);                % Did the subject make the right selection? 1 = yes, 0 = no
            Test_Data.order_of_orientations(k,:) = Preliminary_Data.order_of_orientations(i,:);         % Record the randomized order of image orientations throughout the trial
            Test_Data.log_odds(k) = Preliminary_Data.log_odds(i);
            Test_Data.ratio(k) = Preliminary_Data.ratio(i)
            Test_Data.average_orientations(:,k) = Preliminary_Data.average_orientations(:,i);
            
            Test_Data.current_trial = k;
            
            
            Test_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
            Test_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
            Test_Data.image_length_x = 5;  % Size of the image along x-axis
            Test_Data.image_length_y = 5;
            
            Test_Data.image_template1(k,:) = Preliminary_Data.image_template1(i,:);
            Test_Data.image_template2(k,:) = Preliminary_Data.image_template2(i,:);
            Test_Data.image_template_difference(k,:) = Preliminary_Data.image_template_difference(i,:);
            
            Test_Data.eye_tracker_points{k} = Preliminary_Data.eye_tracker_points{i};
            
            res = Test_Data.screen_resolution;
            Test_Data.left_template = zeros(res * 5);
            for t=1:5
                Test_Data.left_template((t-1)*res+1:t*res, (t-1)*res+1:t*res) = 1;
            end
            Test_Data.right_template = rot90(Test_Data.left_template);
            
            
            image_collection_test(k,:,:,:) = image_collection(i,:,:,:);
                
        k=k+1;
        end
    end
end
if phase == 0       % save the test contrast data
    %% Save final data to folder
    if ~exist(directory, 'dir') % Check the directory actually exists
        mkdir(directory);
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestContrast.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data','image_collection_test'); % save the data
    else
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestContrast.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data','image_collection_test'); % save the data
    end
elseif phase == 1       % save the test ratio data
    %% Save final data to folder
    if ~exist(directory, 'dir') % Check the directory actually exists
        mkdir(directory);
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestRatio.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data','image_collection_test'); % save the data
    else
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborTestRatio.mat']); % create a name for the data you want to save
        save(fileName, 'Test_Data','image_collection_test'); % save the data
    end
end
end