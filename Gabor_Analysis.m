function [] = Gabor_Analysis(subjectID, groupings, preliminary, phase, directory)

%% Analyze All Data
if preliminary == 1 || preliminary == 2
    %load the preliminary data
    if phase == 0
        filename = fullfile(directory, 'RawData', [subjectID '-GaborDataContrast.mat']);
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
        xlabel('Trial'), ylabel('Contrast Level')
        title('Contrast Level by Trial Number')
        
        Get_Figure('Visual Preliminary Reaction Time Analysis');
        %plot reaction time by trial number
        x = 1:Preliminary_Data.current_trial;
        y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
        e = plot(x,y);
        %set(e,'Linewidth',2);
        xlabel('Trial'), ylabel('Reaction Time in msecs')
        title('Reaction Time by Trial Number')
        
    elseif phase == 1
        filename = fullfile(directory, 'RawData', [subjectID '-GaborDataRatio.mat']);
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
        xlabel('Trial'), ylabel('Ratio Level')
        title('Ratio Level by Trial Number')
        
        Get_Figure('Visual Preliminary Reaction Time Analysis');
        %plot reaction time by trial number
        x = 1:Preliminary_Data.current_trial;
        y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
        e = plot(x,y);
        %set(e,'Linewidth',2);
        xlabel('Trial'), ylabel('Reaction Time in msecs')
        title('Reaction Time by Trial Number')
    end
    
    [unique_contrast_conditions, ~, IC] = unique(Preliminary_Data.contrast);
    num_trials_at_x = zeros(length(unique_contrast_conditions),1);
    for i=1:length(unique_contrast_conditions)
        num_trials_at_x(i) = sum(IC == i);
    end
    
    % How many different types trials did the subject see?
    temp_orientation=Preliminary_Data.order_of_orientations';
    temp1=sum(temp_orientation(:,:));
    temp2=Preliminary_Data.number_of_images - temp1;
    temp= temp1-temp2;
    [unique_ratio_conditions,~,IC] = unique(temp);
    
    num_trials_at_x_ratio = zeros(length(unique_ratio_conditions),1);
    for i=1:length(unique_ratio_conditions)
        num_trials_at_x_ratio(i) = sum(IC == i);
    end
    
    average_over_contrast_trials = zeros(1,length(unique_contrast_conditions)); % This will be figuring out how often the subject chose left for each type of trial
    average_over_ratio_trials = zeros(1,length(unique_ratio_conditions)); % This will be figuring out how often the subject chose left for each type of trial
    
    choice = Preliminary_Data.choice;
    accuracy = Preliminary_Data.accuracy;
    
    contrast = Preliminary_Data.contrast;
    % We need to know the contrast values and how many different kind of trials did the subject encounter?
    
    
    contrast_proportions = zeros(1,length(unique_contrast_conditions));
    for i = 1:length(contrast_proportions)
        contrast_proportions(i) = sum(contrast == unique_contrast_conditions(i)); % Check which unique condition does the contrast match
    end
    
    ratio_proportions = zeros(1,length(unique_ratio_conditions));
    for i = 1:length(ratio_proportions)
        ratio_proportions(i) = sum(temp == unique_ratio_conditions(i)); % Check which unique condition does the flash/click rate match
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
        answer = choice(temp == unique_ratio_conditions(p));   % If the trial has the right flash/click rate, get the subject choice
        ratio_stderror(p) = std(answer)/sqrt(length(answer));  % Take the standard error of the choices made
        ratio_yes(p) = sum(choice .* (temp == unique_ratio_conditions(p)));
        average_over_ratio_trials(p) = sum(choice .* (temp == unique_ratio_conditions(p)))/ratio_proportions(p);  % Calculate the average probability of the subject choosing left
    end
    
    
    
    
end

%% Analyze Test Data
if preliminary == 0 || preliminary == 2
    
    
    
    if phase == 0
        % Load data for test gabor experiment
        
        filename = fullfile(directory, 'RawData', [subjectID '-GaborTestContrast.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            return;
        else
            results=load(filename); % Will load Test_Data
        end
        
        
    elseif phase == 1
        % Load data for test gabor experiment
        
        filename = fullfile(directory, 'RawData', [subjectID '-GaborTestRatio.mat']);
        if ~exist(filename, 'file')
            disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
            return;
        else
            results=load(filename);  % Will load Test_Data
        end
        
    end
    %% Serial Dependencies
    [prob_correct_left, prob_correct_right, prob_wrong_left, prob_wrong_right] = Serial_Dependencies(results.Test_Data);     % Print out the serial dependencies
    f=Get_Figure('Serial Dependences: Probability of choosing Left');
    axis off;
    t=uitable(f,'Position',[1 180 500 70],'RowName',{'Correct in prev trial';'Incorrect in prev trial'});
    d={prob_correct_left, prob_correct_right;prob_wrong_left, prob_wrong_right};
    t.Data=d;
    t.ColumnName={'Chose Left in prev trial','Chose Right in prev trial'};
    
    %% Graph Psychometric Curve
    %{
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
    
    %}
    
    addpath(genpath([directory 'Code']));
    Get_Figure('Psychometric Curve And Analysis');
    subplot(2,4,[1,2,5,6]); hold on;  % Plot the psychometric function of the high and low contrast/volume condition
    if phase == 0
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
        
        plotOptions.xLabel         = 'Contrast Level';     % xLabel
        plotOptions.yLabel         = 'Percent Correct';    % yLabel
        legend('Contrast Trials', 'Location','southeast')
        
        
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
        f = errorbar(unique_contrast_conditions,average_over_contrast_trials,contrast_stderror,'bs');  % Black line
        %set(f,'Linewidth',2); 
        hold on;   % Make the lines thicker
        axis tight
        
    elseif phase == 1
        
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
        f = errorbar(unique_ratio_conditions,average_over_ratio_trials,ratio_stderror,'bs');  % Black line
        %set(f,'Linewidth',2); 
        hold on;   % Make the lines thicker
        axis tight
        legend('Ratio Trials', 'Location','southeast')
        
        x = unique_contrast_conditions;
        fitValuesr = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);

        
    end
    
    test_collection_of_images = results.image_collection_test;
    test_image_template_difference = results.Test_Data.image_template_difference;
    test_choice = results.Test_Data.choice;
    
    [trials, number_of_images, ht, wd] = size(test_collection_of_images);
    % number of trials, images shown per trial, and the height and width of the image in the experiment
    
    sublength = number_of_images / groupings;
    
    res=results.Test_Data.screen_resolution;
    h=ht/res;
    w=wd/res;
    
    %========================================
    
    % Image * (T1 - T2)
    im=zeros(trials, h, w, number_of_images);
    for i=1:trials
        for j=1:number_of_images
            v=mat2cell(squeeze(test_collection_of_images(i,j,:,:)), res*ones(1,h), res*ones(1,w));
            for k=1:length(v)
                for l=1:length(v)
                    im(i,k,l,j)=sum(v{k,l}(:))/(res*res);
                end
            end
        end
        
    end
    
    cell_images = mat2cell(im, ones(trials, 1), h, w, number_of_images);
    cell_images = cellfun(@squeeze, cell_images, 'UniformOutput', false);
    [weights, ~, errors] = ...
        CustomRegression.PsychophysicalKernelImage(cell_images, test_choice, 0, 0, 10, 0, 0, 0);




    %X = test_image_template_difference;  % trials x frames
    %Y = test_choice; % 1 x trials
    %[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, Y', [2 number_of_images/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
    % Last term of wAR1 and SDebars are bias terms

    %Get_Figure('Gabor Logistic');
    %subplot(1,1,1); hold on;
    subplot(2,4,[3,4,7,8]);hold on;
    errorbar(1:number_of_images, weights(h*w+1:end-1), errors(h*w+1:end-1));  % Blue plot
    plot([mean(reshape(1:number_of_images, [sublength groupings]))], [sum(reshape(weights(h*w+1:end-1), [sublength groupings]))],'r*-');    % Red plot

    %x=squeeze(mean(reshape(permute(X, [2 1]), [sublength groupings length(choice)])))';
    %[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(x, Y', [2 groupings/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
    %errorbar([0.5:(groupings-0.5)]*sublength, wAR1(1:end), SDebars(1:end),'g-');  % Green plot
    %axis tight;
    %xlabel('Image Frame'), ylabel('Weight')
    title('Weighting the Image Frames')
    %legend('Actual Weight in each frame','Actual weights summed in groups','Weights obtained after averaging groups of images','Location','southoutside')
    legend('Actual Weight in each frame','Actual weights summed in groups','Location','southoutside')




    Get_Figure('Image Kernel');
    r=reshape(weights(1:h*w), [h w]);
    imagesc(r);
    colorbar;



    %}


    if phase ==0
        Get_Figure('Psychometric Curve of ratio for volume/contrast trials');

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
        f = errorbar(unique_ratio_conditions,average_over_ratio_trials,ratio_stderror,'bs');  % Black line
        %set(f,'Linewidth',2); 
        hold on;   % Make the lines thicker
        axis tight
        legend('Ratio Trials', 'Location','southeast')
        
        x = unique_ratio_conditions;
        fitValuesr = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);


    end

    % Plot the absolute values of delta clicks/volume
    Get_Figure('Psychometric Curve of absolute ratio');
    %[params,ll] = fit_sigmoid(flash_rate,Test_Data.choice);
    %p = compute_sigmoid(unique_ratio_conditions,params(1),params(2),params(3));
    unique_ratio_conditions = sort(unique_ratio_conditions);
    x=max(unique_ratio_conditions(unique_ratio_conditions<=0));
    y=min(unique_ratio_conditions(unique_ratio_conditions>=0));
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

    set(best_fit1,'Linewidth',2); hold on;   % Make the lines thicker
    set(best_fit2,'Linewidth',2); hold on;   % Make the lines thicker
    axis([0 inf 0 inf])

    if phase==0
        title('Psychometric Curve of absolute ratio');
        xlabel('|Ratio Difference|')
        ylabel('Percentage of chosing left')
        legend('Absolute Ratio PS for volume trials(<=0)', 'Absolute Ratio PS for contrast trials(>=0)','Location','southoutside')
    end
    if phase==1
        title('Psychometric Curve of absolute ratio');
        xlabel('|Ratio Difference|')
        ylabel('Percentage of chosing left')
        legend('Absolute Ratio PS for ratio trials(<=0)', 'Absolute Ratio PS for ratio trials(>=0)','Location','southoutside')
    end



end

end