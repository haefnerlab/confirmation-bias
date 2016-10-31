function [] = IdealAuditory(subjectID, automatic, phase, add_noise, directory)

% Example Input - Experiment('Matthew', 'Matthew', 0, 0, '/Users/bcs206/Documents/Summer/')

if ~exist('automatic','var')  || ~exist('add_noise','var') || ~exist('directory','var')
    automatic = 0;       % 0 = normal subject, 1 = ideal observer auto running the experiment, 2 = human_like like observer, 3 = human observer
    phase = 0;           % 0 = volume experiment, 1 = ratio experiment
    add_noise = 0;       % 0 = no noise added, 1 = soome background noise added to the sound
    directory = pwd;     % Make it equal to the current directory
end

% This functions initializes and runs an experiment using PsychToolBox

% subjectID is a string to dictate which subject is currently running the
% experiment. Ex. Experiment('01')

% automatic determines if it's to be a person or the computer running the experiment
% Ex. Experiment('01', 0) = person, Experiment('01', 1) = computer auto-runs as the ideal observer

% phase determines if the subject is being tested on volume/contrast or
% ratio. 0 = volume/contrast, 1 = ratio


% add_noise includes some background noise
% Ex. Experiment('01', 0, 0) = no noise added, Experiment('01', 0, 1) = some noise added

% directory allows this code to be able to create and save files of the subject data on any computer


%% Set Up the Initialization of the expeirment
cd([directory 'Code/']) % Set the current directory


    
    fileName = sprintf('%s%s-AuditoryDataVolume.mat',[directory 'RawData/'],subjectID); % Set the desired filename of the experimental data
    if ~exist(fileName, 'file') % Check to see if the subject has already done the preliminary phase or not
        
        preliminary_trials = 10000;
        loops = 1;
        
        Preliminary_Data.move_on = zeros(1,preliminary_trials*loops);          % Is the subject ready to move on or not? Always 0 or 1 for how many trials they got right so far
        Preliminary_Data.step_size = zeros(1,preliminary_trials*loops);        % By how much to adjust the contrast [1.5, 1.2, or 1.1]
        Preliminary_Data.reversal_counter = zeros(1,preliminary_trials*loops);   % How many trials has the subject got wrong? When to change the step size?
        Preliminary_Data.volume = zeros(1,preliminary_trials*loops);         % How loud the sound is, or the signal level
        Preliminary_Data.correct_answer = zeros(1,preliminary_trials*loops);         % What was the right ear/answer?
        Preliminary_Data.reaction_time = zeros(1,preliminary_trials*loops);          % How quick is the subject to input their answer for the orientation choice? Recorded in ms
        Preliminary_Data.choice = zeros(1,preliminary_trials*loops);                 % Which ear did the subject choose? 1 = left, 0 = right
        Preliminary_Data.accuracy = zeros(1,preliminary_trials*loops);               % Did the subject make the right selection? 1 = yes, 0 = no
        Preliminary_Data.ratio = zeros(1,preliminary_trials*loops);  % Re ords underlying click rate
        Preliminary_Data.staircase_answer = zeros(1,preliminary_trials*loops);%stores answers based on actual click rates for staircase
        
        
        Preliminary_Data.current_trial = 0;
        
        % Set properties of a click
        Preliminary_Data.sampling_rate = 6000;   % Indirectly controls frequency
        Preliminary_Data.stimulus_duration = 2;           % How many seconds is each trial?
        Preliminary_Data.bins = Preliminary_Data.sampling_rate * Preliminary_Data.stimulus_duration;         % There needs to be an number of bins for the sound vector
        
        Preliminary_Data.click_rate = zeros(2,preliminary_trials*loops);          % How many clicks did each ear actually hear?
        Preliminary_Data.average_clicks = zeros(2,preliminary_trials*loops);      % Average number of clicks to be the mean of the poisson distribution
        Preliminary_Data.number_of_frames = 60 * Preliminary_Data.stimulus_duration; % How many 'frames' or chances for a click per second in a trial?
        Preliminary_Data.counts = zeros(2,Preliminary_Data.number_of_frames);
        
        Preliminary_Data.order_of_clicks = zeros(preliminary_trials*loops, 2, Preliminary_Data.number_of_frames);
        % Store all clicks sounded, ever
        % Number of trials and the sampling rate times stimulus duration to store the two sounds for the two ears
        
        max_volume = 0.01;      % Starting background volume level
        volume = max_volume;
        
        % Begin Preliminary Trials
        i=1;
        while i <= preliminary_trials * loops
            
            disp(i)
            
            Preliminary_Data.current_trial = i;
            
            
            Preliminary_Data.volume(i) = volume;
            
            if rand < 0.5   % The left ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = 20;
                Preliminary_Data.average_clicks(2,i) = 4;
            else            % The right ear will hear more clicks
                Preliminary_Data.average_clicks(1,i) = 4;
                Preliminary_Data.average_clicks(2,i) = 20;
            end
            
            Preliminary_Data.ratio(i) = 20;
            
            % The correct anxwer is based on underlying click rate, not on the actual number of clicks each ear hears
            wPtr=1;
            I = trialStimuliAuditory(wPtr, subjectID, Preliminary_Data, i, automatic, phase, add_noise, directory);
            
            Preliminary_Data.reaction_time(i) = I.reaction;
            Preliminary_Data.choice(i) = I.choice;     % If 1, subject chose left, and if 0, the subject chose right
            Preliminary_Data.click_rate(1,i) = I.number_of_left_clicks;
            Preliminary_Data.click_rate(2,i) = I.number_of_right_clicks;   % Actual # of clicks for each frame
            
            Preliminary_Data.order_of_clicks(i, :, :) = I.clicks; % Record binary vector indicating the prescence or lack of click in each time bin
            Preliminary_Data.counts (1,:)= Preliminary_Data.counts (1,:)+I.l;
            Preliminary_Data.counts (2,:)= Preliminary_Data.counts (2,:)+I.r;
            
            
            
            %% Save the data after every trial in case of shut-downs
            if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
                mkdir([directory 'RawData/']);
                fileName = sprintf('%s%s-AuditoryDataVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            else
                fileName = sprintf('%s%s-AuditoryDataVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
                save(fileName, 'Preliminary_Data'); % save the data
            end
            i=i+1;
        end
        
        %% Save final data to folder
        if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
            mkdir([directory 'RawData/']);
            fileName = sprintf('%s%s-AuditoryDataVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        else
            fileName = sprintf('%s%s-AuditoryDataVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
            save(fileName, 'Preliminary_Data'); % save the data
        end
        
        %Test_Data=Preliminary_Data;
    end
    %{
    %% Save final data to folder
    if ~exist([directory 'RawData/'], 'dir') % Check the directory actually exists
        mkdir([directory 'RawData/']);
        fileName = sprintf('%s%s-AuditoryTestVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    else
        fileName = sprintf('%s%s-AuditoryTestVolume.mat',[directory 'RawData/'],subjectID); % create a name for the data you want to save
        save(fileName, 'Test_Data'); % save the data
    end
    %}
end   
