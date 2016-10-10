function Experiment_Analysis(opts)
% EXPERIMENT_ANALYSIS Analyze psychometric and psychophysical kernel for a
% particular subject by creating a series of figures.
%
%   EXPERIMENT_ANALYSIS(opts) takes a single struct argument 'opts' with
%   the following fields. The 'subjectID' field is required.
%
%   opts.subjectID [string] - the subject to retrieve the data from. Ex.
%   '01', 'Matthew', etc, etc,...
%
%   opts.automatic [0-1] - if the data  being analyzed will be auditory or
%   visual. Ex. 0 = visual bar task, 1 = auditory click task
%
%   opts.preliminary [0-1] - if we want to skip the analysis of the
%   preliminary data.
%
%   opts.phase [0-2] - if we are analyzing the volume or ratio phase. 0 =
%   skip the preliminary analysis, 1 = skip the test analysis, 2 = don't
%   skip either analysis phase
%
%   opts.manual [0-1] - if we are analyzing manually saved data or not. Ex.
%   0 = automatically saved test data, 1 = manually saved test data
%
%   opts.median [0-2] - if we want to analyse test data above median 0,
%   below median 1 or analyse the entire test data 2.
%
%   opts.version [0-1] - if the weights of the Ideal observer (or manually
%   chosen weights) will be plotted on the graphs along with the weights of
%   the data beign analyzed. Ex. 0 = no ideal weights plotted, 1 = yes plot
%   the ideal weights
%
%   opts.difference [0-1] - how we graph the difference in the weights for
%   the left and right side. 0 means we directly subtract the right weights
%   from the left weights. 1 means we take the difference in the
%   flash/click rates and feed that difference into the regression model
%   and graph that output with errorbars
%
%   opts.directory [string path] - allows this code to be able to create
%   and save files of the subject data on any computer

% Set defaults

if ~isfield(opts, 'subjectID'),   error('opts.subjectID is required'); end
if ~isfield(opts, 'automatic'),   opts.automatic =   0; end
if ~isfield(opts, 'preliminary'), opts.preliminary = 0; end
if ~isfield(opts, 'phase'),       opts.phase =       0; end
if ~isfield(opts, 'manual'),      opts.manual =      0; end
if ~isfield(opts, 'median'),      opts.median =      0; end
if ~isfield(opts, 'version'),     opts.version =     0; end
if ~isfield(opts, 'difference'),  opts.difference =  0; end
if ~isfield(opts, 'directory'),   opts.directory =   '../'; end

%% Analyze Visual Data
if opts.automatic == 0
    if opts.preliminary == 1 || opts.preliminary == 2
        %load the preliminary visual data
        if opts.phase == 0
            prelimFile = [opts.directory 'RawData/' opts.subjectID '-VisualDataContrast.mat'];
            if ~exist(prelimFile, 'file')
                disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
                disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
                return;
            else
                load(prelimFile); % Load Preliminary_Data
            end
            
            Get_Figure('Visual Preliminary Contrast Analysis');
            %plot contrast level by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.contrast(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.contrast(Preliminary_Data.test_phase),'r*.');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);
            xlabel('Trial'), ylabel('Contrast Level')
            title('Contrast Level by Trial Number')
            
            Get_Figure('Visual Preliminary Reaction Time Analysis');
            %plot reaction time by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.reaction_time(Preliminary_Data.test_phase),'r*.');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);
            xlabel('Trial'), ylabel('Reaction Time in msecs')
            title('Reaction Time by Trial Number')
        elseif opts.phase == 1
            % Analyze the ratio data
            prelimFile = [opts.directory 'RawData/' opts.subjectID '-VisualDataRatio.mat'];
            if ~exist(prelimFile, 'file')
                disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
                disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
                return;
            else
                load(prelimFile); % Load Preliminary_Data
            end
            
            Get_Figure('Visual Preliminary Volume Analysis');
            %plot ratio level by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.ratio(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.ratio(Preliminary_Data.test_phase),'r*.');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);
            xlabel('Trial'), ylabel('Ratio Level')
            title('Ratio Level by Trial Number')
            
            Get_Figure('Visual Preliminary Reaction Time Analysis');
            %plot reaction time by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.reaction_time(Preliminary_Data.test_phase),'r*.');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);
            xlabel('Trial'), ylabel('Reaction Time')
            title('Reaction Time by Trial Number')
        end
    end
    
    %% Load data for test visual experiment
    if opts.preliminary == 0 || opts.preliminary == 2
        if opts.phase == 0
            if opts.manual == 0
                filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualTestContrast.mat'];
                if ~exist(filename, 'file')
                    filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualNoisyTest.mat'];  % Are you trying to analyze the older data with the older file name?
                    %else
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            elseif opts.manual == 1
                filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualTestContrastManual.mat'];
                if ~exist(filename, 'file')
                    filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualNoisyTest.mat'];  % Are you trying to analyze the older data with the older file name?
                    %else
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            end
        elseif opts.phase == 1
            if opts.manual == 0
                filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualTestRatio.mat'];
                if ~exist(filename, 'file')
                    filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualNoisyTest.mat'];  % Are you trying to analyze the older data with the older file name?
                    %else
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            elseif opts.manual == 1
                filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualTestRatioManual.mat'];
                if ~exist(filename, 'file')
                    filename = [[opts.directory 'RawData/'] opts.subjectID '-VisualNoisyTest.mat'];  % Are you trying to analyze the older data with the older file name?
                    %else
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            end
        end
    end
    
    %% Analyze Auditory Data
elseif opts.automatic == 1
    if opts.preliminary == 1 || opts.preliminary == 2
        %load the first preliminary auditory data
        if opts.phase == 0
            prelimFile = [opts.directory 'RawData/' opts.subjectID '-AuditoryDataVolume.mat'];
            if ~exist(prelimFile, 'file')
                disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
                disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
                return;
            else
                load(prelimFile); % Load Preliminary_Data
            end
            
            
            %plot volume level by trial number
            Get_Figure('Auditory Preliminary Volume Analysis');
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.volume(1, 1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);hold on;
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.volume(Preliminary_Data.test_phase),'bs');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);hold on;
            xlabel('Trial'), ylabel('Volume Level')
            title('Volume by Trials')
            
            Get_Figure('Auditory Preliminary Reaction Time Analysis');
            %plot reaction time by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);hold on;
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.reaction_time(Preliminary_Data.test_phase),'rs');   % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);hold on;
            xlabel('Trial'), ylabel('Reaction Time in msecs')
            title('Reaction Time by Trial Number')
        elseif opts.phase == 1
            prelimFile = [opts.directory 'RawData/' opts.subjectID '-AuditoryDataRatio.mat'];
            if ~exist(prelimFile, 'file')
                disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
                disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
                return;
            else
                load(prelimFile); % Load Preliminary_Data
            end
            
            Get_Figure('Auditory Preliminary Ratio Analysis');
            %plot ratio level by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.ratio(1, 1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);hold on;
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.ratio(Preliminary_Data.test_phase),'bs');   % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);hold on;
            xlabel('Trial'), ylabel('Ratio Level')
            title('Ratio by Trials')
            
            Get_Figure('Auditory Preliminary Reaction Time Analysis');
            %plot reaction time by trial number
            x = 1:Preliminary_Data.current_trial;
            y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
            e = plot(x,y);
            set(e,'Linewidth',2);hold on;
            p = plot(Preliminary_Data.test_phase, Preliminary_Data.reaction_time(Preliminary_Data.test_phase),'rs');    % Plot the points where we switch from preliminary to test phase
            set(p,'Linewidth',2);hold on;
            xlabel('Trial'), ylabel('Reaction Time in msecs')
            title('Reaction Time by Trial Number')
        end
    end
    
    %% Load data for test auditory experiment
    if opts.preliminary == 0 || opts.preliminary == 2
        if opts.phase == 0
            if opts.manual == 0
                filename = [[opts.directory 'RawData/'] opts.subjectID '-AuditoryTestVolume.mat'];
                if ~exist(filename, 'file')
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    disp(strcat('Maybe the test data is saved under a different name?'));
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            elseif opts.manual == 1
                filename = [[opts.directory 'RawData/'] opts.subjectID '-AuditoryTestVolumeManual.mat'];
                if ~exist(filename, 'file')
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    disp(strcat('Maybe the test data is saved under a different name?'));
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            end
            % Since the analysis for the auditory is identical to the analysis for
            % the visual experiment, I decided to just rename variables rather than
            % write up a very similar code file
            
            Test_Data.contrast = Test_Data.volume;                    % The contrast and volume values will be interchangable
            Test_Data.flash_rate = Test_Data.click_rate;              % How many flashes/clicks was planned per trial
            Test_Data.current_trial = Test_Data.current_trial;        % What was the last trial completed in the experimeent?
            Test_Data.choice = Test_Data.choice;                      % Did the subject chose the left or right bar/ear?
            Test_Data.order_of_flashes = Test_Data.order_of_clicks;   % What was the order of flashes/clicks?
            Test_Data.number_of_images = Test_Data.number_of_frames;  % How many images/bins were there in a trial?
            
            Test_Data.order_of_flashes((isnan(Test_Data.choice)),:,:) = [];
            Test_Data.flash_rate(:,(isnan(Test_Data.choice))) = [];
            Test_Data.contrast(isnan(Test_Data.choice))=[];
            Test_Data.accuracy(isnan(Test_Data.choice)) = [];
            Test_Data.choice(isnan(Test_Data.choice)) = [];
            Test_Data.current_trial = length(Test_Data.choice);
            
            
        elseif opts.phase == 1
            if opts.manual == 0
                filename = [[opts.directory 'RawData/'] opts.subjectID '-AuditoryTestRatio.mat'];
                if ~exist(filename, 'file')
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    disp(strcat('Maybe the test data is saved under a different name?'));
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            elseif opts.manual == 1
                filename = [[opts.directory 'RawData/'] opts.subjectID '-AuditoryTestRatioManual.mat'];
                if ~exist(filename, 'file')
                    disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                    disp(strcat('Maybe the test data is saved under a different name?'));
                    return;
                else
                    load(filename);  % Will load Test_Data
                end
            end
            % Since the analysis for the auditory is identical to the analysis for
            % the visual experiment, I decided to just rename variables rather than
            % write up a very similar code file
            
            Test_Data.contrast = Test_Data.volume;                    % The contrast and volume values will be interchangable
            Test_Data.flash_rate = Test_Data.click_rate;              % How many flashes/clicks actually occurred per trial
            Test_Data.current_trial = Test_Data.current_trial;        % What was the last trial completed in the experimeent?
            Test_Data.choice = Test_Data.choice;                      % Did the subject chose the left or right bar/ear?
            Test_Data.order_of_flashes = Test_Data.order_of_clicks;   % What was the order of flashes/clicks?
            Test_Data.number_of_images = Test_Data.number_of_frames;  % How many images/bins were there in a trial?
            
            Test_Data.order_of_flashes((isnan(Test_Data.choice)),:,:) = [];
            Test_Data.flash_rate(:,(isnan(Test_Data.choice))) = [];
            Test_Data.contrast(isnan(Test_Data.choice))=[];
            Test_Data.accuracy(isnan(Test_Data.choice)) = [];
            Test_Data.choice(isnan(Test_Data.choice)) = [];
            Test_Data.current_trial = length(Test_Data.choice);
            
        end
    end
end

%% Analyze Test Data
if opts.preliminary == 0 || opts.preliminary == 2
    
    m = mean(Test_Data.volume(:));
    Trial_Data = Test_Data;
    if opts.median == 1
        Trial_Data.order_of_flashes(find(Test_Data.contrast(:)> m),:,:) = [];
        Trial_Data.flash_rate(:,find(Test_Data.contrast(:)>m)) = [];
        Trial_Data.contrast(find(Test_Data.contrast(:)>m))=[];
        Trial_Data.accuracy(find(Test_Data.contrast(:)>m)) = [];
        Trial_Data.choice(find(Test_Data.contrast(:)>m)) = [];
        Trial_Data.current_trial = length(find(Test_Data.contrast(:)<=m));
    elseif opts.median == 0
        Trial_Data = Test_Data;
        Trial_Data.order_of_flashes(find(Test_Data.contrast(:)<= m),:,:) = [];
        Trial_Data.flash_rate(:,find(Test_Data.contrast(:)<=m)) = [];
        Trial_Data.contrast(find(Test_Data.contrast(:)<=m))=[];
        Trial_Data.accuracy(find(Test_Data.contrast(:)<=m)) = [];
        Trial_Data.choice(find(Test_Data.contrast(:)<=m)) = [];
        Trial_Data.current_trial = length(find(Test_Data.contrast(:)>m));
    end
    Test_Data = Trial_Data;
    [prob_correct_left, prob_correct_right, prob_wrong_left, prob_wrong_right] = Serial_Dependencies(Test_Data);     % Print out the serial dependencies
    f=Get_Figure('Serial Dependences');
    axis off;
    %t=uitable(f,'Position',[150 180 260 60]);
    t=uitable(f,'Position',[150 180 260 60],'RowName',{'%Correct';'%Incorrect'});
    d={prob_correct_left, prob_correct_right;prob_wrong_left, prob_wrong_right};
    t.Data=d;
    t.ColumnName={'Left','Right'};
    
    
    [unique_contrast_conditions, ~, IC] = unique(Test_Data.contrast);
    num_trials_at_x = zeros(length(unique_contrast_conditions),1);
    for i=1:length(unique_contrast_conditions)
        num_trials_at_x(i) = sum(IC == i);
    end
    %unique_contrast_conditions = sort(unique_contrast_conditions);
    
    % How many different types trials did the subject see?
    
    [unique_ratio_conditions,~,IC] = unique(Test_Data.flash_rate(1,:) - Test_Data.flash_rate(2,:));
    num_trials_at_x_ratio = zeros(length(unique_ratio_conditions),1);
    for i=1:length(unique_ratio_conditions)
        num_trials_at_x_ratio(i) = sum(IC == i);
    end
    %unique_ratio_conditions = sort(unique_ratio_conditions);
    
    % How many different types trials did the subject see?
    
    % For example was there a trial with 20 more flashes/clicks on the left vs the right?
    
    average_over_contrast_trials = zeros(1,length(unique_contrast_conditions)); % This will be figuring out how often the subject chose left for each type of trial
    average_over_ratio_trials = zeros(1,length(unique_ratio_conditions)); % This will be figuring out how often the subject chose left for each type of trial
    
    choice = Test_Data.choice;
    accuracy = Test_Data.accuracy;
    
    contrast = Test_Data.contrast;
    % We need to know the contrast values and how many different kind of trials did the subject encounter?
    
    flash_rate = Test_Data.flash_rate(1,:) - Test_Data.flash_rate(2,:);
    % There's a set number of flashes/clicks for the left bar/ear as opposed to the right bar/ear
    % We need to know the difference in the flash rate and how many different kind of trials did the subject encounter?
    
    % How many trials of each type was there?
    contrast_proportions = zeros(1,length(unique_contrast_conditions));
    for i = 1:length(contrast_proportions)
        contrast_proportions(i) = sum(contrast == unique_contrast_conditions(i)); % Check which unique condition does the contrast match
    end
    
    ratio_proportions = zeros(1,length(unique_ratio_conditions));
    for i = 1:length(ratio_proportions)
        ratio_proportions(i) = sum(flash_rate == unique_ratio_conditions(i)); % Check which unique condition does the flash/click rate match
    end
    
    
    contrast_stderror = zeros(1,length(unique_contrast_conditions));  % After getting the number of trials for each type, we need the standard error
    contrast_yes = zeros(1,length(unique_contrast_conditions));
    contrast_no = zeros(1,length(unique_contrast_conditions));
    for p = 1:length(unique_contrast_conditions)
        answer = accuracy(contrast == unique_contrast_conditions(p));   % If the trial has the right contrast level, get the subject choice
        contrast_stderror(p) = std(answer)/sqrt(length(answer));  % Take the standard error of the choices made
        contrast_yes(p) = sum(answer);
        contrast_no(p) = length(answer)-contrast_yes(p);
        average_over_contrast_trials(p) = sum(accuracy .* (contrast == unique_contrast_conditions(p)))/contrast_proportions(p);  % Calculate the average probability of the subject being correct
    end
    
    ratio_stderror = zeros(1,length(unique_ratio_conditions));  % After getting the number of trials for each type, we need the standard error
    ratio_yes = zeros(1,length(unique_ratio_conditions));
    for p = 1:length(unique_ratio_conditions)
        answer = choice(find(flash_rate == unique_ratio_conditions(p)));   % If the trial has the right flash/click rate, get the subject choice
        ratio_stderror(p) = std(answer)/sqrt(length(answer));  % Take the standard error of the choices made
        ratio_yes(p) = sum(choice .* (flash_rate == unique_ratio_conditions(p)));
        average_over_ratio_trials(p) = sum(choice .* (flash_rate == unique_ratio_conditions(p)))/ratio_proportions(p);  % Calculate the average probability of the subject choosing left
    end
    
    order_of_flashes = [squeeze(Test_Data.order_of_flashes(:,1,:)) squeeze(Test_Data.order_of_flashes(:,2,:))];
    order_of_flashes = reshape(order_of_flashes, Test_Data.current_trial, []);
    % Take the order of flashes/clicks in each trial and reorganize it so each trial has a vector of the left side followed by the right side
    % We needed to change order_of_flashes from a 3d matrix into a 2d matrix for the regression methods
    
    %{
    Get_Figure('Psychometric Curve');
    subplot(2,4,[1,2,5,6]); hold on;  % Plot the psychometric function of the high and low contrast/volume condition
    e = errorbar(unique_conditions, average_over_trials, stderror, '.-','color', [0 0 0]); hold on;  % Black line
    set(e,'Linewidth',2); hold on;   % Make the lines thicker
    axis tight
    title('Psychometric Functions');
    
    if opts.phase == 0
        if opts.automatic == 0
            xlabel('Contrast')
            ylabel('Percentage correct')
            legend('Contrast Trials', 'Location','southeast')
        elseif opts.automatic == 1
            xlabel('Volume')
            ylabel('Percentage correct')
            legend('Volume Trials', 'Location','southeast')
        end
    elseif phase == 1
        xlabel('Ratio Difference')
        ylabel('Percentage of chosing left')
        legend('Ratio Trials', 'Location','southeast')
    end
    %}
    
    
    
    %{
    Get_Figure('Psychometric Curve And Analysis');
    subplot(2,4,[1,2,5,6]); hold on;  % Plot the psychometric function of the high and low contrast/volume condition
   
    if opts.phase == 0
        
        % Plot the volume/contrast level by the percentage correct
        [params,ll] = fit_sigmoid(Test_Data.contrast,accuracy);
        p=compute_sigmoid(unique_contrast_conditions,params(1),params(2),params(3));
        e=plot(unique_contrast_conditions,p,'g-');
        set(e,'Linewidth',2); hold on;
        f = errorbar(unique_contrast_conditions,average_over_contrast_trials,contrast_stderror,'bs');  % Black line
        set(f,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        title('Psychometric Functions');
        if opts.automatic == 0
            xlabel('Contrast')
            ylabel('Percentage Correct')
            legend('Contrast Trials', 'Location','southeast')
        elseif opts.automatic == 1
            xlabel('Volume')
            ylabel('Percentage Correct')
            legend('Volume Trials', 'Location','southeast')
        end
    elseif opts.phase == 1
        [params,ll] = fit_sigmoid(flash_rate,Test_Data.choice);
        p=compute_sigmoid(unique_ratio_conditions,params(1),params(2),params(3));
        e=plot(unique_contrast_conditions,p,'g-');hold on;
        e = errorbar(unique_contrast_conditions,average_over_contrast_trials,contrast_stderror,'bs');  % Black line
        set(e,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        title('Psychometric Functions');
        xlabel('Ratio Difference')
        ylabel('Percentage of chosing left')
        legend('Ratio Trials', 'Location','southeast')
    end
    
    %}
    
    %addpath(genpath('psignifit-master'));
    addpath(genpath([opts.directory 'Code']));
    Get_Figure('Psychometric Curve And Analysis');
    subplot(2,4,[1,2,5,6]); hold on;  % Plot the psychometric function of the high and low contrast/volume condition
    if opts.phase == 0
        data_cont_vol = [unique_contrast_conditions(:), contrast_yes(:), num_trials_at_x(:)];
        options=struct;
        options.sigmoidName    = 'weibull';
        options.expType        = '2AFC';
        options.estimateType = 'MLE';
        options.nblocks = length(unique_contrast_conditions);
        options.confP          = [0.95,0.9,.68];
        options.instantPlot    = 0;
        options.setBordersType = 0;
        options.maxBorderValue = .00001;
        options.moveBorders    = 1;
        options.dynamicGrid    = 0;
        options.widthalpha     = .05;
        options.threshPC       = .7;
        options.CImethod       = 'percentiles';
        options.gridSetType    = 'cumDist';
        options.fixedPars      = nan(5,1);
        options.useGPU         = 0;
        options.poolMaxGap     = inf;
        options.poolMaxLength  = inf;
        options.poolxTol       = 0;
        options.betaPrior      = 10;
        options.verbose        = true;
        options.stimulusRange  = 0;
        options.fastOptim      = false;
        
        result=psignifit(data_cont_vol,options);
        
        title('Psychometric Functions');
        plotOptions.dataColor      = [0,105/255,170/255];  % color of the data
        plotOptions.plotData       = 0;                    % plot the data?
        plotOptions.lineColor      = [0,0,0];              % color of the PF
        plotOptions.lineWidth      = 2;                    % lineWidth of the PF
        if opts.automatic == 0
            plotOptions.xLabel         = 'Contrast Level';     % xLabel
            plotOptions.yLabel         = 'Percent Correct';    % yLabel
            legend('Contrast Trials', 'Location','southeast')
        elseif opts.automatic == 1
            plotOptions.xLabel         = 'Volume Level';     % xLabel
            plotOptions.yLabel         = 'Percent Correct';    % yLabel
            legend('Volume Trials', 'Location','southeast')
            
        end
        
        plotOptions.labelSize      = 15;                   % font size labels
        plotOptions.fontSize       = 10;                   % font size numbers
        plotOptions.fontName       = 'Helvetica';          % font
        plotOptions.tufteAxis      = false;                % use special axis
        plotOptions.plotAsymptote  = false;                 % plot Asympotes
        plotOptions.plotThresh     = false;                 % plot Threshold Mark
        plotOptions.aspectRatio    = false;                % set aspect ratio
        plotOptions.extrapolLength = .2;                   % extrapolation percentage
        plotOptions.CIthresh       = false;                % draw CI on threhold
        plotOptions.dataSize       = 10000./sum(result.data(:,3)); % size of the data-dots
        
        h = plotPsych(result,plotOptions);hold on;
        
        x = unique_contrast_conditions;
        fitValuesc = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);
        
        %e=plot(unique_contrast_conditions,fitValuesc,'g-');
        %set(e,'Linewidth',2); hold on;
        f = errorbar(unique_contrast_conditions,average_over_contrast_trials,contrast_stderror,'bs');  % Black line
        set(f,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        
        
        fileName = sprintf('%s%s-result.mat',[opts.directory 'RawData/'],opts.subjectID); % create a name for the data you want to save
        save(fileName, 'result');
        %plotPsych(result,plotOptions);
        
    elseif opts.phase == 1
        
        data_ratio = [unique_ratio_conditions(:), ratio_yes(:), num_trials_at_x_ratio(:)];
        
        options=struct;
        options.sigmoidName    = 'norm';
        options.expType        = 'YesNo';
        options.estimateType = 'MLE';
        options.nblocks = length(unique_contrast_conditions);
        options.confP          = [0.95,0.9,.68];
        options.instantPlot    = 0;
        options.setBordersType = 0;
        options.maxBorderValue = .00001;
        options.moveBorders    = 1;
        options.dynamicGrid    = 0;
        options.widthalpha     = .05;
        options.threshPC       = .7;
        options.CImethod       = 'percentiles';
        options.gridSetType    = 'cumDist';
        options.fixedPars      = nan(5,1);
        options.useGPU         = 0;
        options.poolMaxGap     = inf;
        options.poolMaxLength  = inf;
        options.poolxTol       = 0;
        options.betaPrior      = 10;
        options.verbose        = true;
        options.stimulusRange  = 0;
        options.fastOptim      = false;
        
        result=psignifit(data_ratio,options);
        
        title('Psychometric Functions');
        plotOptions.dataColor      = [0,105/255,170/255];  % color of the data
        plotOptions.plotData       = 0;                    % plot the data?
        plotOptions.lineColor      = [0,0,0];              % color of the PF
        plotOptions.lineWidth      = 2;                    % lineWidth of the PF
        plotOptions.xLabel         = 'Ratio Difference';     % xLabel
        plotOptions.yLabel         = 'Percent Correct';    % yLabel
        plotOptions.labelSize      = 15;                   % font size labels
        plotOptions.fontSize       = 10;                   % font size numbers
        plotOptions.fontName       = 'Helvetica';          % font
        plotOptions.tufteAxis      = false;                % use special axis
        plotOptions.plotAsymptote  = false;                 % plot Asympotes
        plotOptions.plotThresh     = false;                 % plot Threshold Mark
        plotOptions.aspectRatio    = false;                % set aspect ratio
        plotOptions.extrapolLength = .2;                   % extrapolation percentage
        plotOptions.CIthresh       = false;                % draw CI on threhold
        plotOptions.dataSize       = 10000./sum(result.data(:,3)); % size of the data-dots
        
        h = plotPsych(result,plotOptions);hold on;
        
        
        %title('Psychometric Functions');
        x = unique_ratio_conditions;
        fitValuesr = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);
        
        %e=plot(unique_ratio_conditions,fitValuesr,'g-');
        %set(e,'Linewidth',2); hold on;
        f = errorbar(unique_ratio_conditions,average_over_ratio_trials,ratio_stderror,'bs');  % Black line
        set(f,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        
        
        %xlabel('Ratio Difference')
        %ylabel('Percentage Chose Left')
        legend('Ratio Trials', 'Location','southeast')
        
        
    end
    
    %============
    %% Regularized logistic regression
    %addpath(genpath('logisticRegression'));  % Jake's logistic functions
    
    subplot(2,4,[3, 4]); hold on;    % Plot left and right weights for the high contrast/volume case
    colors='br';
    LRWeightsErrors = cell(2,2);
    for LeftRight=[1 2]
        st_idx = (LeftRight-1)*120+1;
        end_idx = st_idx+120-1;
        X = order_of_flashes(:,st_idx:end_idx);
        X=X/std(X(:));
        X=[X ones(size(X,1),1)];  % Add a bias term
        Y = Test_Data.choice(:); % 1 x trials
        wmle=glmfit(X(:,1:end-1), Y, 'binomial');
        %     [wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
        %         Y, [2 Test_Data.number_of_images], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]',wmle); % Call on Jake's functions
        
        %[wAR1,~,SDebars,~] = autoRegress_logisticRidge(X, ...
        %  Y, [1 Test_Data.number_of_images*2], 0.01, 10.^(0:6)',wmle); % Call on Jake's functions
        
        [wAR1,~,SDebars,~] = autoRegress_logisticAR1(X, Y, Test_Data.number_of_images, 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]',wmle); % Call on Jake's functions
        wAR1(end)=[];      % Remove bias term
        SDebars(end)=[];
        SDebars_diff = SDebars;   % Save this for the difference in weights
        LRWeightsErrors{LeftRight,1} = wAR1;
        LRWeightsErrors{LeftRight,2} = SDebars;
        
        % Plot the subject's left weights
        linestyle = ['.-' colors(LeftRight)];
        e = errorbar(wAR1, SDebars, linestyle);
        set(e,'Linewidth',2); hold on;
    end
    legend('Left', 'Right');
    % To print out the subject's/Ideal Observer's left bar weights, uncomment the following:
    % wAR1(1:Test_Data.number_of_images)
    % SDebars(1:Test_Data.number_of_images)
    
    % This is to include Ideal Weights in the graph as a comparison to the subject's weights
    % Since the Ideal Weights will require a different analysis, the weights will have to be manually included
    
    % Ideal Left Weights
    if opts.version == 1
        idealLeftBarWeights = [0.6492,0.6763,0.5486,0.4267,0.4868,0.5528,0.5854,0.6241,0.6171, ...
            0.6158,0.6609,0.7116,0.7902,0.8748,0.9297,0.9909,0.9033,0.8220, ...
            0.7521,0.6882,0.6627,0.6428,0.6599,0.6830,0.6610,0.6452,0.7221, ...
            0.8053,0.7785,0.7577,0.7805,0.8092,0.7482,0.6933,0.7429,0.7985, ...
            0.7910,0.7894,0.6454,0.5072];
        idealLeftBarErrors = [0.4635,0.2257,0.3517,0.2307,0.3438,0.2094,0.3435,0.2287,0.3511, ...
            0.2313,0.3509,0.2340,0.3508,0.2262,0.3512,0.2333,0.3590,0.2498, ...
            0.3566,0.2306,0.3524,0.2337,0.3501,0.2242,0.3460,0.2258,0.3531, ...
            0.2392,0.3545,0.2354,0.3526,0.2320,0.3540,0.2291,0.3508,0.2491, ...
            0.3626,0.2393,0.3531,0.2416];
        e = errorbar(idealLeftBarWeights, idealLeftBarErrors,'.-c');  % Cyan Line
        set(e,'Linewidth',2); hold on;
    end
    
    % Plot the subject's right weights
% % % %     e = errorbar(wAR1(Test_Data.number_of_images+1:end), SDebars(Test_Data.number_of_images+1:end), '.-r');  % Red Line
% % % %     set(e,'Linewidth',2); hold on;
% % % %     legend('Left weights', 'Right weights','Location','northoutside')
    % To print out the subject's/Ideal Observer's right bar weights, uncomment the following:
    % wAR1(Test_Data.number_of_images+1:end)
    % SDebars(Test_Data.number_of_images+1:end)
    
    % Ideal Right Weights
    if opts.version == 1
        idealRightBarWeights = [-0.5964,-0.6116,-0.7484,-0.8809,-0.8447,-0.8040,-0.7814,-0.7543,-0.7903, ...
            -0.8222,-0.7914,-0.7567,-0.6735,-0.5858,-0.5238,-0.4571,-0.5405,-0.6194, ...
            -0.6902,-0.7566,-0.7869,-0.8131,-0.7870,-0.7566,-0.7593,-0.7573,-0.6906, ...
            -0.6193,-0.6487,-0.6737,-0.6526,-0.6272,-0.6799,-0.7282,-0.6615,-0.5903, ...
            -0.5987,-0.6028,-0.7201,-0.8332];
        idealRightBarErrors = [0.4645,0.2279,0.3549,0.2391,0.3552,0.2392,0.3598,0.2485,0.3588, ...
            0.2434,0.3541,0.2368,0.3516,0.2317,0.3481,0.2225,0.3514,0.2437, ...
            0.3530,0.2312,0.3541,0.2401,0.3530,0.2311,0.3483,0.2274,0.3498, ...
            0.2333,0.3521,0.2282,0.3529,0.2425,0.3652,0.2558,0.3543,0.2336, ...
            0.3547,0.2313,0.3552,0.2558];
        e = errorbar(idealRightBarWeights, idealRightBarErrors, '.-m');  % Magenta Line
        set(e,'Linewidth',2); hold on;
    end
    if opts.phase == 0
        title('Low Volume/Contrast PK');
    elseif opts.phase == 1
        title('High Volume/Contrast PK');
    end
    axis tight
    xlabel('Frames')
    ylabel('Weights')
    
    
    subplot(2,4,[7,8]); hold on; % Graph the difference in weights between the left and right bar/ear
    
    if opts.difference == 0
        % Graph the difference in the weights for the left and right side
        wdiff = LRWeightsErrors{1,1}-LRWeightsErrors{2,1};
        SDebars = LRWeightsErrors{1,2}+LRWeightsErrors{2,2};
        %s = SDebars(1:length(wAR1)/2) + SDebars(length(wAR1)/2+1:end);
        %e = errorbar([1:Test_Data.number_of_images], w', s', '.-k');
        e= plot(wdiff, '.-k');
        set(e,'Linewidth',2); hold on;
    else
        % Take the difference in the flash/click rates and feed the difference into the regresssion model
        % and then graph the output with errorbars
        difference_in_stimuli = squeeze(Test_Data.order_of_flashes(:,1,:)) - squeeze(Test_Data.order_of_flashes(:,2,:));
        X = difference_in_stimuli;
        X = X/std(X(:));
        X = [X ones(size(X,1),1)];  % Add a bias term
        Y = Test_Data.choice; % 1 x trials
        %wmle=glmfit(X(:,1:end-1), Y', 'binomial');
%         [wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
%             Y', [2 Test_Data.number_of_images/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
        
        %[wAR1,~,SDebars,~] = autoRegress_logisticRidge(X, ...
        %   Y', [1 Test_Data.number_of_images], 0.01, 10.^(0:6)',wmle);
        [wAR1,~,SDebars,~] = autoRegress_logisticAR1(X, ...
           Y', Test_Data.number_of_images, 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]',wmle);
        wAR1(end)=[];
        SDebars(end)=[];
        e = errorbar(wAR1(1:end), SDebars(1:end),'.-k');
        set(e,'Linewidth',2); hold on;
    end
    
    title('PK Difference (Left Image Weights - Right Image Weights)');
    axis tight
    xlabel('Frames')
    ylabel('Weights')
    if opts.phase == 0
        legend('Volume','Location','southoutside')
    elseif opts.phase == 1
        legend('Ratio','Location','southoutside')
    end
    
    
    %% Other additional plots
    
    
    % Plot the change in change in click/flash rate by the percentage left for the threshold volume/contrast case
    if opts.phase ==0
        Get_Figure('Psychometric Curve of ratio for volume/contrast trials');
        %{
        h = plot(unique_ratio_conditions, ratio_proportions./sum(ratio_proportions), 'k-');hold on;
        set(h,'Linewidth',2); hold on;   % Make the lines thicker
        
        [params,ll] = fit_sigmoid(flash_rate,Test_Data.choice);
        p = compute_sigmoid(unique_ratio_conditions,params(1),params(2),params(3));
        best_fit = plot(unique_ratio_conditions,p,'g-');hold on;
        set(best_fit,'Linewidth',2); hold on;   % Make the lines thicker
        e = errorbar(unique_ratio_conditions,average_over_ratio_trials,ratio_stderror,'bs'); hold on;  % Black line
        set(e,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        title('Psychometric Function with ratio difference for Volume/Contrast Trials');
        xlabel('Ratio Difference')
        ylabel('Percentage of chosing left')
        %}
        data_ratio = [unique_ratio_conditions(:), ratio_yes(:), num_trials_at_x_ratio(:)];
        
        options=struct;
        options.sigmoidName    = 'norm';
        options.expType        = 'YesNo';
        options.estimateType = 'MLE';
        options.nblocks = length(unique_contrast_conditions);
        options.confP          = [0.95,0.9,.68];
        options.instantPlot    = 0;
        options.setBordersType = 0;
        options.maxBorderValue = .00001;
        options.moveBorders    = 1;
        options.dynamicGrid    = 0;
        options.widthalpha     = .05;
        options.threshPC       = .7;
        options.CImethod       = 'percentiles';
        options.gridSetType    = 'cumDist';
        options.fixedPars      = nan(5,1);
        options.useGPU         = 0;
        options.poolMaxGap     = inf;
        options.poolMaxLength  = inf;
        options.poolxTol       = 0;
        options.betaPrior      = 10;
        options.verbose        = true;
        options.stimulusRange  = 0;
        options.fastOptim      = false;
        
        result=psignifit(data_ratio,options);
        
        title('Psychometric Functions');
        plotOptions.dataColor      = [0,105/255,170/255];  % color of the data
        plotOptions.plotData       = 0;                    % plot the data?
        plotOptions.lineColor      = [0,0,0];              % color of the PF
        plotOptions.lineWidth      = 2;                    % lineWidth of the PF
        plotOptions.xLabel         = 'Ratio Difference in Volume trials';     % xLabel
        plotOptions.yLabel         = 'Percent Correct';    % yLabel
        plotOptions.labelSize      = 15;                   % font size labels
        plotOptions.fontSize       = 10;                   % font size numbers
        plotOptions.fontName       = 'Helvetica';          % font
        plotOptions.tufteAxis      = false;                % use special axis
        plotOptions.plotAsymptote  = false;                 % plot Asympotes
        plotOptions.plotThresh     = false;                 % plot Threshold Mark
        plotOptions.aspectRatio    = false;                % set aspect ratio
        plotOptions.extrapolLength = .2;                   % extrapolation percentage
        plotOptions.CIthresh       = false;                % draw CI on threhold
        plotOptions.dataSize       = 10000./sum(result.data(:,3)); % size of the data-dots
        
        h = plotPsych(result,plotOptions);hold on;
        
        
        %title('Psychometric Functions');
        x = unique_ratio_conditions;
        fitValuesr = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);
        
        %e=plot(unique_ratio_conditions,fitValuesr,'g-');
        %set(e,'Linewidth',2); hold on;
        f = errorbar(unique_ratio_conditions,average_over_ratio_trials,ratio_stderror,'bs');  % Black line
        set(f,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        
        
        %xlabel('Ratio Difference')
        %ylabel('Percentage Chose Left')
        legend('Ratio Trials', 'Location','southeast')
        
        
        
    end
    
    % Plot the absolute values of delta clicks/volume
    Get_Figure('Psychometric Curve of absolute ratio');
    %[params,ll] = fit_sigmoid(flash_rate,Test_Data.choice);
    %p = compute_sigmoid(unique_ratio_conditions,params(1),params(2),params(3));
    unique_ratio_conditions = sort(unique_ratio_conditions);
    x=max(unique_ratio_conditions(find(unique_ratio_conditions<=0)));
    y=min(unique_ratio_conditions(find(unique_ratio_conditions>=0)));
    for i=1:length(unique_ratio_conditions)
        if unique_ratio_conditions(i)==x
            inflection_point1 = i;
        end
        if unique_ratio_conditions(i)==y
            inflection_point2 = i;
        end
    end
    
    u1=abs(unique_ratio_conditions(1:inflection_point2-1));
    u1(end+1)=-unique_ratio_conditions(inflection_point2);
    p1=fitValuesr(1:inflection_point2-1);
    p1(end+1)=fitValuesr(inflection_point2);
    best_fit1 = plot(u1,p1,'b-');
    best_fit2 = plot(unique_ratio_conditions(inflection_point1:end),fitValuesr(inflection_point1:end),'r-');
    size(fitValuesr)
    size(unique_ratio_conditions)
    set(best_fit1,'Linewidth',2); hold on;   % Make the lines thicker
    set(best_fit2,'Linewidth',2); hold on;   % Make the lines thicker
    axis([0 inf 0 inf])
    
    if opts.phase==0
        title('Psychometric Curve of absolute ratio');
        xlabel('|Ratio Difference|')
        ylabel('Percentage of chosing left')
        legend('Absolute Ratio PS for volume trials(<=0)', 'Absolute Ratio PS for volume trials(>=0)','Location','southoutside')
    end
    if opts.phase==1
        title('Psychometric Curve of absolute ratio');
        xlabel('|Ratio Difference|')
        ylabel('Percentage of chosing left')
        legend('Absolute Ratio PS for ratio trials(<=0)', 'Absolute Ratio PS for ratio trials(>=0)','Location','southoutside')
    end
    
end


end

