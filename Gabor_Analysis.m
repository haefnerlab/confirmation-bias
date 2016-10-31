function [] = Gabor_Analysis(subjectID_prelim,subjectID, groupings, preliminary, phase, manual, directory)

%% Analyze All Data
if preliminary == 1 || preliminary == 2
    %load the preliminary data
    if phase == 0
        filename = fullfile(directory, 'RawData', [subjectID_prelim '-GaborDataContrast.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
            return;
        else
            load(filename); % Load Preliminary_Data
        end
        
        Get_Figure('Preliminary Analysis');
        
        %plot contrast level by trial number
        x = 1:Preliminary_Data.current_trial;
        y = Preliminary_Data.contrast(1:Preliminary_Data.current_trial);
        plot(x,y)
        plot(Preliminary_Data.test_phase, Preliminary_Data.contrast(Preliminary_Data.test_phase),'rs');    % Plot the points where we switch from preliminary to test phase
        xlabel('Trial'), ylabel('Contrast Level')
        title('Contrast Level by Trial Number')
        
    elseif phase == 1
        filename = fullfile(directory, 'RawData', [subjectID_prelim '-GaborDataRatio.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
            return;
        else
            load(filename); % Load Preliminary_Data
        end
        
        Get_Figure('Preliminary Analysis');
        
        %plot ratio level by trial number
        x = 1:Preliminary_Data.current_trial;
        y = Preliminary_Data.ratio(1:Preliminary_Data.current_trial);
        plot(x,y);
        plot(Preliminary_Data.test_phase, Preliminary_Data.ratio(Preliminary_Data.test_phase),'bs');    % Plot the points where we switch from preliminary to test phase
        xlabel('Trial'), ylabel('Ratio Level')
        title('Ratio Level by Trial Number')
    end
end

%% Analyze Test Data
if preliminary == 0 || preliminary == 2
    
    Serial_Dependencies(Test_Data);     % Print out the serial dependencies
    
    if phase == 0
        % Load data for test gabor experiment
        if manual == 0
            filename = fullfile(directory, 'RawData', [subjectID '-GaborTestContrast.mat']);
            if ~exist(filename, 'file')
                disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                return;
            else
                results = load(filename); % Will load Test_Data
            end
        elseif manual == 1
            filename = fullfile(directory, 'RawData', [subjectID '-GaborTestContrastManual.mat']);
            if ~exist(filename, 'file')
                disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                return;
            else
                results = load(filename); % Will load Test_Data
            end
        end
    elseif phase == 1
        % Load data for test gabor experiment
        if manual == 0
            filename = fullfile(directory, 'RawData', [subjectID '-GaborTestRatio.mat']);
            if ~exist(filename, 'file')
                disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                return;
            else
                results = load(filename);  % Will load Test_Data
            end
        elseif manual == 1
            filename = fullfile(directory, 'RawData', [subjectID '-GaborTestRatioManual.mat']);
            if ~exist(filename, 'file')
                disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
                return;
            else
                results = load(filename);  % Will load Test_Data
            end
        end
    end
    
    %% Graph Psychometric Curve
    for i=1:results.Test_Data.current_trial
        orientation_rate(1,i) = sum(results.Test_Data.order_of_orientations(i,:));
        orientation_rate(2,i) = results.Test_Data.number_of_images - orientation_rate(1,i);
    end
    if phase == 0
        unique_conditions = unique(results.Test_Data.contrast);
        % How many different types trials did the subject see?
    elseif phase == 1
        unique_conditions = unique(orientation_rate(1,:) - orientation_rate(2,:));
        % How many different types trials did the subject see?
        % For example was there a trial with 5 more left orientation images vs the right orientation images?
    end
    
    average_over_trials = zeros(1,length(unique_conditions)); % This will be figuring out how often the subject chose left for each type of trial
    
    choice = results.Test_Data.choice;
    accuracy = results.Test_Data.accuracy;
    
    if phase == 0
        contrast = results.Test_Data.contrast;
        % We need to know the contrast values and how many different kind of trials did the subject encounter?
    elseif phase == 1
        orientation_diff = orientation_rate(1,:) - orientation_rate(2,:);
        % There's a set number of orientations for the left image as opposed to the right image
        % We need to know the difference in the orientation rate and how many different kind of trials did the subject encounter?
    end
    
    % How many trials of each type was there?
    proportions = zeros(1,length(unique_conditions));
    for i = 1:length(proportions)
        if phase == 0
            proportions(i) = sum(contrast == unique_conditions(i)); % Check which unique condition does the contrast match
        elseif phase == 1
            proportions(i) = sum(orientation_diff == unique_conditions(i)); % Check which unique condition does the orientation rate match
        end
    end
    
    stderror = zeros(1,length(unique_conditions));  % After getting the number of trials for each type, we need the standard error
    
    for p = 1:length(unique_conditions)
        if phase == 0
            answer = accuracy(find(contrast == unique_conditions(p)));   % If the trial has the right contrast level, get the subject choice
            stderror(p) = std(answer)/sqrt(length(answer));  % Take the standard error of the choices made
            average_over_trials(p) = sum(accuracy .* (contrast == unique_conditions(p)))/proportions(p);  % Calculate the average probability of the subject choosing left
        elseif phase == 1
            answer = choice(find(orientation_diff == unique_conditions(p)));   % If the trial has the right orientation rate, get the subject choice
            stderror(p) = std(answer)/sqrt(length(answer));  % Take the standard error of the choices made
            average_over_trials(p) = sum(choice .* (orientation_diff == unique_conditions(p)))/proportions(p);  % Calculate the average probability of the subject choosing left
        end
    end
    %{
    x = unique_conditions;
    y = average_over_trials;
    sigm_function = @(A, x)(A(1) ./ (1 + exp(-A(2)*x)));
    A0 = ones(1,2);
    [A_fit,~,~,~,~,ErrorModelInfo] = nlinfit(x, y, sigm_function, A0, 'ErrorModel','proportional','ErrorParameters',0.5);
    
    Get_Figure('Best Fit Psychometric Curves');
    subplot(2,4,[1,2,5,6]); hold on;  % Plot the psychometric function of the high and low contrast/volume condition
    % Best Fit
    e = errorbar(x, (A_fit(1) ./ (1 + exp(-A_fit(2)*x))), ErrorModelInfo.ErrorParameters, '.-','color', [0 0 0]); hold on;  % Black line
    set(e,'Linewidth',2); hold on;   % Make the lines thicker
    % Data Points
    e = plot(unique_conditions, average_over_trials, 'b*'); hold on;  % Blue line
    set(e,'Linewidth',2); hold on;   % Make the lines thicker
    axis tight
    title('Best Fit Psychometric Function');
    ylabel('Percentage of Choosing Left')
    %}
    Get_Figure('Psychometric Curve And Gabor Logistic');
    subplot(2,4,[1,2,5,6]); hold on; 
    choice = results.Test_Data.choice;
    accuracy = results.Test_Data.accuracy;
    contrast = results.Test_Data.contrast;
    if phase == 0
        [params,ll] = fit_sigmoid(contrast,accuracy);
        p=compute_sigmoid(unique_conditions,params(1),params(2),params(3));
        %e=plot(unique_conditions,p);
        e = plot(unique_conditions,average_over_trials,'bs',unique_conditions,p,'g-');hold on;  % Black line
        set(e,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        title('Psychometric Functions');
        
        xlabel('Contrast')
        ylabel('Percentage correct')
        legend('Contrast Trials', 'Location','southeast')
        
    elseif phase == 1
        [params,ll] = fit_sigmoid(orientation_diff,choice);
        p=compute_sigmoid(unique_conditions,params(1),params(2),params(3));
        %e=plot(unique_conditions,p);
        e = plot(unique_conditions,average_over_trials,'bs',unique_conditions,p,'g-');hold on;  % Black line
        set(e,'Linewidth',2); hold on;   % Make the lines thicker
        axis tight
        title('Psychometric Functions');
        xlabel('Ratio Difference')
        ylabel('Percentage of chosing left')
        legend('Ratio Trials', 'Location','southeast')
    end
    
    
    
    
    
    
    
    
    collection_of_images = results.test_image_collection;
    image_template_difference = results.Test_Data.image_template_difference;
    choice = results.Test_Data.choice;
    
    [~, number_of_images, ~, ~] = size(collection_of_images);
    % number of trials, images shown per trial, and the height and width of the image in the experiment
    
    sublength = number_of_images / groupings;
    
    
   
    
    %========================================
    
    % Image * (T1 - T2)
    
    
    X = image_template_difference;  % trials x frames
    Y = choice; % 1 x trials
    [wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, Y', [2 number_of_images/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
    % Last term of wAR1 and SDebars are bias terms
    
    %Get_Figure('Gabor Logistic');
    %subplot(1,1,1); hold on;
    subplot(2,4,[3,4,7,8]);hold on;
    errorbar([1:number_of_images], wAR1(1:end), SDebars(1:end));  % Blue plot
    plot([mean(reshape(1:number_of_images, [sublength groupings]))], [sum(reshape(wAR1(1:end), [sublength groupings]))],'r*-');    % Red plot
    
    x=squeeze(mean(reshape(permute(X, [2 1]), [sublength groupings length(choice)])))';
    [wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(x, Y', [2 groupings/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
    errorbar([0.5:(groupings-0.5)]*sublength, wAR1(1:end), SDebars(1:end),'g-');  % Green plot
    axis tight;
    xlabel('Image Frame'), ylabel('Weight')
    title('Weighting the Image Frames')
    legend('Actual Weight in each frame','Actual weights summed in groups','Weights obtained after averaging groups of images','Location','southoutside')
    
end
end