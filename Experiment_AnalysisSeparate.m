function [] = Experiment_AnalysisSeparate(subjectID_prelim, subjectID, automatic, version, difference, average, directory)

% Example Input - Experiment_AnalysisSeparate('Matthew', 0, 0, 0, 0, 0, '/Users/bcs206/Documents/Summer/')


% This function is an analysis of the visual/auditory data produced by the
% Bar Task and the Poisson Click Task

% subjectID determines the subject to retrieve the data from. Ex. '01', 'Matthew', etc, etc,...

% automatic determines if the data  being analyzed will be auditory or
% visual. Ex. 0 = visual bar task, 1 = auditory click task

% version determines if the weights of the Ideal observer (or manually chosen weights) will be plotted
% on the graphs along with the weights of the data beign analyzed.
% Ex. 0 = no ideal weights plotted, 1 = yes plot the ideal weights

% difference dictates how we graph the difference in the weights for the left and right side
% 0 means we directly subtract the right weights from the left weights
% 1 means we take the difference in the flash/click rates and feed that difference into the regression model and graph that output with errorbars

% average dictates how to take the average of the difference in the weights
% 0 means we get and plot 4 averages of the difference in the weights before the difference is fed into the regression methods
% 1 means we get and plot 4 averages of the difference in the weights after the difference is fed into the regression methods

% directory allows this code to be able to create and save files of the subject data on any computer


if ~exist('automatic','var') || ~exist('version','var')  % Check for missing arguments
    automatic = 0;
    version = 0;
	directory = '/Users/bcs206/Documents/Summer/';
end

if automatic == 0
	%load the preliminary visual data
    prelimFile = fullfile(directory, 'RawData', [subjectID_prelim '-VisualPreliminary.mat']);
    if ~exist(prelimFile, 'file')
		prelimFile = fullfile(directory, 'RawData', [subjectID_prelim '-VisualNoisyPreliminary.mat']); % Are you trying to analyze the older data with the older file name?
		if ~exist(prelimFile, 'file')
			disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
			disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
			return;
		else	
			load(prelimFile); % Load Preliminary_Data
		end
    else
		load(prelimFile); % Load Preliminary_Data
	end
    
    Get_Figure('Visual Preliminary Analysis A');
         %plot step size by trial number
    x = 1:Preliminary_Data.current_trial;
    y = Preliminary_Data.step_size(1:Preliminary_Data.current_trial);
    plot(x,y)
    xlabel('Trial'), ylabel('Step Size')
    title('Step Size by Trial Number')
    
	Get_Figure('Visual Preliminary Analysis B');
        %plot contrast level by trial number
    x = 1:Preliminary_Data.current_trial;
    y = Preliminary_Data.contrast(1:Preliminary_Data.current_trial);
    plot(x,y)
    xlabel('Trial'), ylabel('Contrast Level')
    title('Contrast Level by Trial Number')
    
	Get_Figure('Visual Preliminary Analysis C');
        %plot reaction time1 by trial number
    x = 1:Preliminary_Data.current_trial;
    y = Preliminary_Data.reaction_time(1:Preliminary_Data.current_trial);
    plot(x,y)
    xlabel('Trial'), ylabel('RT')
    title('Reaction Time by Trial Number')

	
    %% Load data for test visual experiment
    
    filename = fullfile(directory, 'RawData', [subjectID '-VisualTest.mat']);
    if ~exist(filename, 'file')
		filename = fullfile(directory, 'RawData', [subjectID '-VisualNoisyTest.mat']);  % Are you trying to analyze the older data with the older file name?
	%else
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        return;
    else
        load(filename);  % Will load Test_Data
    end
    	
elseif automatic == 1
	%load the first preliminary auditory data
    prelimFile = fullfile(directory, 'RawData', [subjectID_prelim '-AuditoryPreliminaryVolume.mat']);
    if ~exist(prelimFile, 'file')
		prelimFile = fullfile(directory, 'RawData', [subjectID_prelim '-AuditoryNoisyPreliminary.mat']); % Are you trying to analyze the older data with the older file name?
		if ~exist(prelimFile, 'file')
			disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
			disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
		else
			load(prelimFile); % Load Preliminary_Data
		end
    else
		load(prelimFile); % Load Preliminary_Data
	end
	
    Get_Figure('Auditory Preliminary Volume Analysis');
		%plot volume level by trial number
	x = 1:Preliminary_Data.current_trial;
	y = Preliminary_Data.volume(1, 1:Preliminary_Data.current_trial);
	e = plot(x,y);
	set(e,'Linewidth',2);
	xlabel('Trial'), ylabel('Volume Level')
	title('Volume by Trials')
    
    
    %load the second preliminary auditory data
    prelimFile = fullfile(directory, 'RawData', [subjectID_prelim '-AuditoryPreliminaryRatio.mat']);
    if ~exist(prelimFile, 'file')
        disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
        disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
    else
        load(prelimFile); % Load Preliminary_Data
    end
	
    Get_Figure('Auditory Preliminary Ratio Analysis');
		%plot ratio by trial number
	x = 1:Preliminary_Data.current_trial;
	y = Preliminary_Data.ratio(1, 1:Preliminary_Data.current_trial);
	e = plot(x,y);
	set(e,'Linewidth',2);
	xlabel('Trial'), ylabel('Underlying Click Rate')
	title('Underlying Click Rate by Trials')
	
    
    %% Load data for test auditory experiment
    filename = [fullfile(directory, 'RawData'] subjectID '-AuditoryTest.mat');
    if ~exist(filename, 'file')
		filename = [fullfile(directory, 'RawData'] subjectID '-AuditoryNoisyTest.mat'); % Are you trying to analyze the older data with the older file name?
	%else
        disp(strcat('ERROR! Missing File: ', filename));  % Return an error message for missing file
        return;
    else
        load(filename);  % Will load Test_Data
    end
    
    % Since the analysis for the auditory is identical to the analysis for
    % the visual experiment, I decided to just rename variables rather than
    % write up a very similar code file
    
    Test_Data.flash_rate = Test_Data.click_rate;              % How many flashes/clicks was planned per trial
    Test_Data.current_trial = Test_Data.current_trial;        % What was the last trial completed in the experimeent?
    Test_Data.choice = Test_Data.choice;                      % Did the subject chose the left or right bar/ear?
    Test_Data.order_of_flashes = Test_Data.order_of_clicks;   % What was the order of flashes/clicks?
    Test_Data.number_of_images = Test_Data.number_of_frames;  % How many images/bins were there in a trial?
end

unique_conditions = unique(Test_Data.flash_rate(1,:) - Test_Data.flash_rate(2,:));
		% How many different types trials did the subject see?
		% For example was there a trial with 20 more flashes/clicks on the left vs the right?

average_over_white_trials = zeros(1,length(unique_conditions)); % This will be figuring out how often the subject chose left for each type of trial
average_over_gray_trials = zeros(1,length(unique_conditions));
	% White refers to the high contrast/volume case
	% Gray refers to the low contrast/volume case

highContrastIndex = logical([ones(1, Test_Data.current_trial/2) zeros(1, Test_Data.current_trial/2)]');
lowContrastIndex = ~highContrastIndex;
	% It's assumed that the high contrast/volume case is always the first half of the trials in the test phase

white_choice = Test_Data.choice(highContrastIndex);  % Get the subject choice for the high and low contrast/volume trials
gray_choice = Test_Data.choice(lowContrastIndex);

white_flash_rate = Test_Data.flash_rate(1,highContrastIndex) - Test_Data.flash_rate(2,highContrastIndex);
gray_flash_rate = Test_Data.flash_rate(1,lowContrastIndex) - Test_Data.flash_rate(2, lowContrastIndex);
	% There's a set number of flashes/clicks for the left bar/ear as opposed to the right bar/ear
	% We need to know the difference in the flash rate and how many different kind of trials did the subject encounter?

% How many trials of each type was there?
white_proportions = zeros(1,length(unique_conditions));
for i = 1:length(white_proportions)
    white_proportions(i) = sum(white_flash_rate == unique_conditions(i)); % Check which unique condition does the flash/click rate match
end
gray_proportions = zeros(1,length(unique_conditions));
for i = 1:length(gray_proportions)
    gray_proportions(i) = sum(gray_flash_rate == unique_conditions(i));
end

white_stderror = zeros(1,length(unique_conditions));  % After getting the number of trials for each type, we need the standard error
gray_stderror = zeros(1,length(unique_conditions));

for p = 1:length(unique_conditions)
    answer = zeros(1,white_proportions(p));   % answer will be a vector of 1s and 0s to correspond to the subject choice in a specific subset of trials
    for temp = 1:Test_Data.current_trial/2
        if white_flash_rate(temp) == unique_conditions(p)
            answer(temp) = white_choice(temp);   % If the trial has the right flash/click rate, get the subject choice
        end
    end
    white_stderror(p) = std(answer);  % Take the standard error of the choices made
    
    answer = zeros(1,gray_proportions(p));
    for temp = 1:Test_Data.current_trial/2
        if gray_flash_rate(temp) == unique_conditions(p)
            answer(temp) = gray_choice(temp);
        end
    end
    gray_stderror(p) = std(answer);
    
    average_over_white_trials(p) = sum(white_choice .* (white_flash_rate == unique_conditions(p)))./white_proportions(p);  % Calculate the average probability of the subject choosing left
    average_over_gray_trials(p) = sum(gray_choice .* (gray_flash_rate == unique_conditions(p)))./gray_proportions(p);
end
order_of_flashes = [squeeze(Test_Data.order_of_flashes(:,1,:)) squeeze(Test_Data.order_of_flashes(:,2,:))];
order_of_flashes = reshape(order_of_flashes, Test_Data.current_trial, []);
	% Take the order of flashes/clicks in each trial and reorganize it so each trial has a vector of the left side followed by the right side


Get_Figure('High Contrast Bar Plot');
	% Plot the psychometric function of the high contrast/volume condition
e = errorbar(unique_conditions, average_over_white_trials, white_stderror, '.-','color', [0 0 0]); hold on;  % Black line
set(e,'Linewidth',2); hold on;   % Make the lines thicker
axis tight
title('High Contrast Psychometric Function');
xlabel('Probability Distribution')
ylabel('Percentage of Choosing Left')

Get_Figure('Low Contrast Bar Plot');
	% Plot the psychometric function of the high contrast/volume condition
e = errorbar(unique_conditions, average_over_gray_trials, gray_stderror, '.-','color', [0.5 0.5 0.5]);  % Gray line
set(e,'Linewidth',2); hold on;
axis tight
title('Low Contrast Psychometric Function');
xlabel('Probability Distribution')
ylabel('Percentage of Choosing Left')

%============
% Regularized logistic regression
addpath('logisticRegression/regressionTools')  % Jake's logistic functions

Get_Figure('RegressionA');
	% Plot left and right weights for the high contrast/volume case
ix=highContrastIndex;
X = order_of_flashes(ix,:);
X=X/std(X(:));
X=[X ones(size(X,1),1)];  % Add a bias term
Y = Test_Data.choice(ix); % 1 x trials
[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
Y', [2 Test_Data.number_of_images], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]'); % Call on Jake's functions
wAR1(end)=[];      % Remove bias term
SDebars(end)=[];
wAR1_white_diff = wAR1;   % Save this for the difference in weights

% Plot the subject's left weights
e = errorbar(wAR1(1:Test_Data.number_of_images), SDebars(1:Test_Data.number_of_images),'.-b'); % Blue Line
set(e,'Linewidth',2); hold on;

% This is to include Ideal Weights in the graph as a comparison to the subject's weights
% Since the Ideal Weights will require a different analysis, the weights will have to be manually included

% Ideal Left Weights
if version == 1
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
e = errorbar(wAR1(Test_Data.number_of_images+1:end), SDebars(Test_Data.number_of_images+1:end), '.-r');  % Red Line
set(e,'Linewidth',2); hold on;

% Ideal Right Weights
if version == 1
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
title('High Amplitude PK');
axis tight
xlabel('Time in Frames')
ylabel('Weights')

Get_Figure('RegressionB');
	% Plot the weights for the low contrast/volume case
ix=lowContrastIndex;
X = order_of_flashes(ix,:);
X=X/std(X(:));
X=[X ones(size(X,1),1)];  % Add a bias term
Y = Test_Data.choice(ix); % 1 x trials
[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
Y', [2 Test_Data.number_of_images], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]'); % Call on Jake's logistic function
wAR1(end)=[];     % Remove bias term
SDebars(end)=[];
wAR1_gray_diff = wAR1;  % Save this for the difference in weights graph

% Plot subject's left weights
e = errorbar(wAR1(1:Test_Data.number_of_images), SDebars(1:Test_Data.number_of_images),'.-b');
set(e,'Linewidth',2); hold on;

% Plot Ideal Left Weights
if version == 1
    idealLeftBarWeights = [1.3359,1.5579,1.3009,1.2588,1.1026,1.1612,1.1775,1.4087,1.2541, ...
        1.3144,1.1093,1.1187,0.9614,1.0186,0.8143,0.8245,0.9267,1.2434, ...
        1.2050,1.3813,1.1789,1.1912,1.0550,1.1338,1.0078,1.0963,1.2744, ...
        1.6670,1.3347,1.2172,1.0856,1.1686,1.0471,1.1406,1.1921,1.4588, ...
        1.2039,1.1641,0.8971,0.8447];
    idealLeftBarErrors = [0.9525,0.4031,0.7042,0.4149,0.6914,0.3764,0.7058,0.4177,0.7094, ...
        0.3847,0.7005,0.3850,0.6982,0.3950,0.6936,0.3540,0.7102,0.4432, ...
        0.7205,0.3764,0.6994,0.4054,0.7169,0.4469,0.7340,0.4748,0.7301, ...
        0.4384,0.7140,0.4203,0.7019,0.4000,0.7156,0.4644,0.7390,0.4630, ...
        0.7192,0.3717,0.6984,0.3952];
    e = errorbar(idealLeftBarWeights, idealLeftBarErrors,'.-c');
    set(e,'Linewidth',2); hold on;
end

% Plot subject's right weights
e = errorbar(wAR1(Test_Data.number_of_images+1:end), SDebars(Test_Data.number_of_images+1:end),'.-r');
set(e,'Linewidth',2); hold on;

% Plot Ideal right weights
if version == 1
    idealRightBarWeights = [-0.9495,-0.7395,-1.0086,-1.0627,-1.2278,-1.1782,-1.1743,-0.9555,-1.1112, ...
        -1.0522,-1.2626,-1.2586,-1.4163,-1.3596,-1.5493,-1.5245,-1.4140,-1.0891, ...
        -1.1324,-0.9611,-1.1646,-1.1534,-1.3005,-1.2328,-1.3628,-1.2784,-1.0931, ...
        -0.6932,-1.0207,-1.1333,-1.2687,-1.1894,-1.3054,-1.2062,-1.1465,-0.8716, ...
        -1.1082,-1.1298,-1.3664,-1.3885];
    idealRightBarErrors = [0.9529,0.4119,0.7043,0.4159,0.6983,0.4062,0.7107,0.4140,0.7239, ...
        0.4443,0.7237,0.4139,0.7210,0.4505,0.7254,0.4136,0.7119,0.3988, ...
        0.7079,0.3795,0.6912,0.3726,0.7049,0.4376,0.7148,0.4157,0.6916, ...
        0.3642,0.6871,0.3938,0.6909,0.3929,0.7149,0.4681,0.7262,0.4149, ...
        0.7166,0.4171,0.7240,0.4401];
    e = errorbar(idealRightBarWeights, idealRightBarErrors, '.-m');
    set(e,'Linewidth',2); hold on;
end
title('Low Amplitude PK');
axis tight
xlabel('Time in Frames')
ylabel('Weights')

% Use differenct legends depending on whether or not we are including Ideal weights in the graph
if version == 1
    legend('Subject Left Bar Weights','Ideal Left Bar Weights','Subject Right Bar Weights','Ideal Right Bar Weights','Location','southoutside')
else
    legend('Left Bar Weights','Right Bar Weights','Location','southoutside')
end

Get_Figure('DifferenceA');
	% Graph the difference in weights between the left and right bar/ear

if difference == 0
	% Graph the difference in the weights for the left and right side
	wAR1 = wAR1_white_diff;
	white_w = wAR1(1:length(wAR1)/2) - wAR1(length(wAR1)/2+1:end);
	e = plot([1:Test_Data.number_of_images], white_w','.-k');
	set(e,'Linewidth',2); hold on;
else
	% Take the difference in the flash/click rates and feed the difference into the regresssion model
	% and then graph the output with errorbars
	difference_in_stimuli = squeeze(Test_Data.order_of_flashes(:,1,:)) - squeeze(Test_Data.order_of_flashes(:,2,:));
	X = difference_in_stimuli(highContrastIndex,:);
	X = X/std(X(:));
	X = [X ones(size(X,1),1)];  % Add a bias term
	Y = Test_Data.choice(highContrastIndex); % 1 x trials
	[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
	Y', [2 Test_Data.number_of_images/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
	wAR1(end)=[];
	white_w = wAR1;
	SDebars(end)=[];
	e = errorbar(wAR1(1:end), SDebars(1:end),'.-k');
	set(e,'Linewidth',2); hold on;
end
title('PK High Contrast Difference (Left Image Weights - Right Image Weights)');
axis tight
xlabel('Time in Frames')
ylabel('Weights')

Get_Figure('DifferenceB');
	% Graph the difference in weights between the left and right bar/ear

if difference == 0
	% Graph the difference in the weights for the left and right side
	wAR1 = wAR1_gray_diff;
	gray_w = wAR1(1:length(wAR1)/2) - wAR1(length(wAR1)/2+1:end);
	e = plot([1:Test_Data.number_of_images], gray_w','.-','Color',[0.5 0.5 0.5]);
	set(e,'Linewidth',2); hold on;
else
	% Take the difference in the flash/click rates and feed the difference into the regresssion model
	% and then graph the output with errorbars
	difference_in_stimuli = squeeze(Test_Data.order_of_flashes(:,1,:)) - squeeze(Test_Data.order_of_flashes(:,2,:));
	X = difference_in_stimuli(lowContrastIndex,:);
	X = X/std(X(:));
	X = [X ones(size(X,1),1)];  % Add a bias term
	Y = Test_Data.choice(lowContrastIndex); % 1 x trials
	[wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
	Y', [2 Test_Data.number_of_images/2], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]');
	wAR1(end)=[];
	gray_w = wAR1;
	SDebars(end)=[];
	e = errorbar(wAR1(1:end), SDebars(1:end),'.-','Color',[0.5 0.5 0.5]);
	set(e,'Linewidth',2); hold on;
end
title('PK Low Contrast Difference (Left Image Weights - Right Image Weights)');
axis tight
xlabel('Time in Frames')
ylabel('Weights')


Get_Figure('AveragesA');
	% Completely new figure to graph the averages
if average == 0;  % We are finding the average of the difference in the weights before it's fed into the regression methods
	plot([mean(reshape(1:length(white_w), [length(white_w)/4 4]))], [sum(reshape(white_w(1:end), [length(white_w)/4 4]))], '.-','color', [0 0 0]);
	% white_w is a variable I'm reusing from the above graph so I don't need to repeat the calculations for the difference in weights
else
	% We are finding the average of the difference in the weights after it's fed into the regression methods
	plot([mean(reshape(1:length(white_w), [length(white_w)/4 4]))], [sum(reshape(white_w(1:end), [length(white_w)/4 4]))], '.-','color', [0 0 0]);
	% white_wAR1 is a variable I'm reusing from the above graph so I don't need to repeat the calculations for the difference in weights after the regression
end
title('Average High Contrast Weight Difference');
axis tight
xlabel('Time in Frames')
ylabel('Weights')


Get_Figure('AveragesB');
	% Completely new figure to graph the averages
if average == 0;  % We are finding the average of the difference in the weights before it's fed into the regression methods
	plot([mean(reshape(1:length(gray_w), [length(gray_w)/4 4]))], [sum(reshape(gray_w(1:end), [length(gray_w)/4 4]))], '.-','color', [0.5 0.5 0.5]);
	% gray_w is a variable I'm reusing from the above graph so I don't need to repeat the calculations for the difference in weights
else
	% We are finding the average of the difference in the weights after it's fed into the regression methods
	plot([mean(reshape(1:length(gray_w), [length(gray_w)/4 4]))], [sum(reshape(gray_w(1:end), [length(gray_w)/4 4]))], '.-','color', [0.5 0.5 0.5]);
	% gray_wAR1 is a variable I'm reusing from the above graph so I don't need to repeat the calculations for the difference in weights after the regression
end
title('Average Low Contrast Weight Difference');
axis tight
xlabel('Time in Frames')
ylabel('Weights')

end