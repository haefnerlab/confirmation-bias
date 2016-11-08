function [] = IdealGabor(subjectID, directory)


cd(fullfile(directory, 'Code')) % Set the current directory



fileName = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']); % Set the desired filename of the experimental data
if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
    
    preliminary_trials = 1600;
    loops = 1;
    
    Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
    Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
    Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
    Preliminary_Data.contrast = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
    Preliminary_Data.number_of_images = 10;
    Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
    Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
    Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
    Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
    Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
    Preliminary_Data.order_of_orientations = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);        % Record the randomized order of image orientations throughout the trial
    Preliminary_Data.log_odds = zeros(1,preliminary_trials*loops);
    Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);
    Preliminary_Data.average_orientations = zeros(2,preliminary_trials*loops);
    
    Preliminary_Data.current_trial = 0;
    
    
    Preliminary_Data.screen_frame = 12;	% how long each image will show on screen in frame rates
    Preliminary_Data.screen_resolution = 25;          % how many pixels correspond to a single datapoint of a gabor
    Preliminary_Data.image_length_x = 5;  % Size of the image along x-axis
    Preliminary_Data.image_length_y = 5;
    
    Preliminary_Data.image_template1 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
    Preliminary_Data.image_template2 = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
    Preliminary_Data.image_template_difference = zeros(preliminary_trials*loops,Preliminary_Data.number_of_images);
    
    
    
    Preliminary_Data.left_template = eye(Preliminary_Data.image_length_y, Preliminary_Data.image_length_x);
    Preliminary_Data.right_template = rot90(Preliminary_Data.left_template);
    
    
    image_collection = zeros(preliminary_trials*loops, Preliminary_Data.number_of_images, ...
        Preliminary_Data.image_length_x, Preliminary_Data.image_length_y);
    
    contrast = 64;
    ratio_sum = Preliminary_Data.number_of_images;
    ratio = 2;
    
    % Begin Preliminary Trials
    i=1;
    while i <= preliminary_trials * loops
        
        disp(i)
        
        Preliminary_Data.current_trial = i;
        
        
        Preliminary_Data.contrast(i) = contrast;
        if rand < 0.5   % The left ear will hear more clicks
            Preliminary_Data.average_orientations(1,i) = ratio;
            Preliminary_Data.average_orientations(2,i) = ratio_sum - ratio;
            desired_orientations = [1, 0];        % Left Orientation
            probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];		% Calculates % of images for one orientation vs another orientation
            Preliminary_Data.correct_answer(i) = 1;	% Since left appears >= than the right orientation, it's the correct answer
            Preliminary_Data.ratio(i) = ratio;  % Use only the one orientation
        else            % The right ear will hear more clicks
            Preliminary_Data.average_orientations(1,i) = ratio_sum - ratio;
            Preliminary_Data.average_orientations(2,i) = ratio;
            desired_orientations = [0, 1];        % Right Orientation
            probabilities = [ratio/ratio_sum, (ratio_sum - ratio)/ratio_sum];		% Calculates % of images for one orientation vs another orientation
            Preliminary_Data.correct_answer(i) = 0;	% Since right appears >= than the left orientation, it's the correct answer
            Preliminary_Data.ratio(i) = ratio;  % Use only the one orientation
        end
        
        
        [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, Preliminary_Data.number_of_images);
        Preliminary_Data.correct_answer(i) = correct_orientation_answer;
        Preliminary_Data.choice(i) = Preliminary_Data.correct_answer(i);
        image_array = makeImages(Preliminary_Data, order_of_orientations, contrast);
        
        % Store all images shown
        image_collection(i,:,:,:) = image_array;
        
        Preliminary_Data.order_of_orientations(i,:) = order_of_orientations;  % Record random ordering of all orientations
        
    i=i+1;    
    end
    if ~exist(directory, 'dir') % Check the directory actually exists
        mkdir(directory);
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']); % create a name for the data you want to save
        save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
    else
        fileName = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']); % create a name for the data you want to save
        save(fileName, 'Preliminary_Data', 'image_collection'); % save the data
    end
end

end