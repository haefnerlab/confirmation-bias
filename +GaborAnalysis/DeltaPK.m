function [perSubjectFigs, combinedFig] = DeltaPK(subjectIDs, phases, per_subject_plots, method, datadir)
%GABORANALYSIS.DELTAPK creates one figure per subject and a combined figure
%(if 2 or more subjects) showing temporal psychophysical kernel analysis.
%
% GABORANALYSIS.DELTAPK(subjectIDs, phases, datadir)
%
% Inputs:
% - subjectIDs: a cell array of subject IDs (strings)
% - phases: 0 for contrast, 1 for ratio, 2 for noise. If 2 numbers given,
%           the difference phases(1)-phases(2) is plotted.
% - datadir: (optional) override the default place to look for data files.

if nargin < 3, per_subject_plots = false; end
if nargin < 5, method = 'reg-lr'; end
if nargin < 6, datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

window_low = 0.5;
window_high = 0.7;

color1 = [32 74 135] / 255;
color2 = [164 0 0] / 255;
color3 = [0 0 0] / 255;
white = [255 255 255] / 255;

% Confidence interval
conf = 0.68;

%% Helper functions

    function trials = getSubjectThresholdTrials(subjectId, phase)
        stair_var = get_stair_var(phase);
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        trials = SubjectData.(stair_var) <= thresh & SubjectData.(stair_var) >= floor;
    end

    function true_pk = getGroundTruthModelPK(subjectId, phase)
        stair_var = get_stair_var(phase);
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        
        if isfield(SubjectData, 'model_pk')
            true_pk = SubjectData.model_pk;
        else
            true_pk = nan(1, SubjectData.number_of_images);
        end
    end

    function [boot_params, frames] = getSubjectBootstrapRegression(subjectId, phase, model)
        stair_var = get_stair_var(phase);
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        frames = SubjectData.number_of_images;
        warning off;
        [floor, thresh] = GaborAnalysis.getThresholdWindow(subjectId, phase, window_low, window_high, datadir);
        warning on;
        % Hard reset of thresholds in the 'ratio' task since bins are too coarse. Performance is
        % within a few pecent of 70% for all subjects at the 4:6 or 6:4 ratios.
        if phase == 1
            floor = .4;
            thresh = .6;
        end
        SubjectDataThresh = GaborThresholdTrials(SubjectData, phase, thresh, floor);
        switch lower(model)
            case {'lr-reg', 'reg-lr'}
                % Regularized logistic regression
                
                % UNCOMMENT TO RUN CROSS-VALIDATION PER SUBJECT
                % memo_name = ['PK-xValid-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                % nFold = 10;
                % hprs = [0 logspace(-3, 5, 9)];
                % [hprs, ~] = LoadOrRun(@CustomRegression.xValidatePK, ...
                %     {SubjectDataThresh.ideal_frame_signals, SubjectDataThresh.choice == +1, hprs, 0, hprs, 1, nFold}, ...
                %     fullfile(memodir, memo_name));
                
                % USING EMPIRICAL GEOMETRIC MEAN OF CROSS-VALIDATED HYPERPARAMETERS
                if phase == 1
                    hprs = [0.1 0 10];
                elseif phase  == 2
                    hprs = [0 0 500];
                end
                
                % Get logistic regression weights with bootstrapping
                memo_name = ['Boot-PK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '-reg.mat'];
                [~, ~, ~, ~, boot_params] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, 500, hprs, 0, true}, ...
                    fullfile(memodir, memo_name));
            case {'regress', 'logistic', 'lr'}
                % Unregularized logistic regression
                hprs = [0 0 0];
                memo_name = ['Boot-PK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, boot_params] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, 500, hprs, 0, true}, ...
                    fullfile(memodir, memo_name));
            case {'cta'}
                % Choice-triggered average
                memo_name = ['Boot-CTA-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, boot_params] = LoadOrRun(@BootstrapCTA, ...
                    {SubjectDataThresh, 500}, ...
                    fullfile(memodir, memo_name));
            case {'lin', 'linear'}
                % Linear model: weight = slope * frameno + offset
                memo_name = ['Boot-lin-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, ~, boot_params] = LoadOrRun(@BootstrapLinearPKFit, ...
                    {SubjectDataThresh, 500, 0, true}, ...
                    fullfile(memodir, memo_name));
            case {'exp', 'exponential'}
                % Exponential model: weight = alpha * exp(beta * frameno)
                memo_name = ['Boot-exp-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, ~, boot_params] = LoadOrRun(@BootstrapExponentialWeightsGabor, ...
                    {SubjectDataThresh, 500, 0, true}, ...
                    fullfile(memodir, memo_name));
            otherwise
                error('Unknown PK method: %s', model);
        end
    end

    function [weights, bias] = paramsToKernel(params, modelstring, frames, normalize)
        switch lower(modelstring)
            case {'regress', 'logistic', 'lr', 'reg-lr', 'lr-reg'}
                weights = params(:,1:end-1);
                bias = params(:,end);
            case {'cta'}
                weights = params;
                bias = zeros(size(params,1), 1);
            case {'lin', 'linear'}
                % In linear model, 3 parameters are [slope offset bias]
                weights = (0:frames-1).*params(:,1) + params(:,2);
                bias = params(:,3);
            case {'exp', 'exponential'}
                % In exponential model, 3 parameters are [alpha beta bias]
                weights = params(:,1).*exp(params(:,2).*(0:frames-1));
                bias = params(:,3);
            otherwise
                error('Unrecognized PK model type: %s', modelstring);
        end
        
        if normalize
            weights = weights ./ mean(weights, 2);
            % NOTE we don't touch 'bias' here!
        end
    end

%% Loop over subjects. Compute individual and combined kernels.

combinedFig = -1;
if length(subjectIDs) > 1
    combinedFig = figure;
    hold on;
end

% CombinedKernelsByPhase{i} contains all bootstrapped temporal weights, normalized to mean 1, per
% phase per subject. They will later be concatenated to get the average across all subjects.
PerSubjectKernelsByPhase = cell(length(phases), length(subjectIDs));

perSubjectFigs = [];
for i=1:length(subjectIDs)
    if per_subject_plots
        perSubjectFigs(i) = figure;
        hold on;
    end
    subjectId = subjectIDs{i};
    if length(phases) == 1
        [boot_params, frames] = getSubjectBootstrapRegression(subjectId, phases, method);
        [boot_weights, boot_bias] = paramsToKernel(boot_params, method, frames, true);
        [~, lo_w, hi_w, med_w] = meanci(boot_weights, conf);
        [~, lo_b, hi_b, med_b] = meanci(boot_bias, conf);
        
        PerSubjectKernelsByPhase{i} = boot_weights;
        
        % Add to the combined figure
        if length(subjectIDs) > 1
            figure(combinedFig);
            plot(1:frames, med_w, 'LineWidth', 1, 'Color', (color1+white)/2);
        end
        
        if per_subject_plots
            boundedline(1:frames, med_w', [med_w-lo_w; hi_w-med_w]');
            errorbar(frames+1, med_b, med_b-lo_b, hi_b-med_b, 'LineWidth', 2, 'Color', color1);
            trials = getSubjectThresholdTrials(subjectId, phases);
            title([strrep(get_stair_var(phases), '_', ' ') 'temporal kernel (' num2str(sum(trials)) '/' num2str(length(trials)) ')']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 frames], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [0 1]);
            xlabel('Time (s)');
            ylabel('Psychophysical Kernel');
            axis tight;
            
            ground_truth_pk = getGroundTruthModelPK(subjectId, phase);
            if ~all(isnan(ground_truth_pk))
                plot(1:frames, ground_truth_pk, 'Color', color1, 'LineWidth', 2);
            end
        end
    elseif length(phases) == 2
        [boot_params1, frames1] = getSubjectBootstrapRegression(subjectId, phases(1), method);
        [boot_w1, boot_b1] = paramsToKernel(boot_params1, method, frames1, true);
        [~, lo_w1, hi_w1, med_w1] = meanci(boot_w1, conf);
        [~, lo_b1, hi_b1, med_b1] = meanci(boot_b1, conf);
        PerSubjectKernelsByPhase{1, i} = boot_w1;
        
        [boot_params2, frames2] = getSubjectBootstrapRegression(subjectId, phases(2), method);
        [boot_w2, boot_b2] = paramsToKernel(boot_params2, method, frames2, true);
        [~, lo_w2, hi_w2, med_w2] = meanci(boot_w2, conf);
        [~, lo_b2, hi_b2, med_b2] = meanci(boot_b2, conf);
        PerSubjectKernelsByPhase{2, i} = boot_w2;

        assert(frames1 == frames2, 'Cannot subtract PKs with different # frames');
        boot_diff_w = boot_w1 - boot_w2;
        [~, lo_diff, hi_diff, med_diff] = meanci(boot_diff_w, conf);
        
        % Add to the combined figure
        % Combined chematic is 3 subplots: each pk in their own and their difference in the third
        if length(subjectIDs) > 1
            figure(combinedFig);
            subplot(1,3,1); hold on;
            plot(1:frames1, med_w1, 'LineWidth', 1, 'Color', (color1+white)/2);
            
            subplot(1,3,2); hold on;
            plot(1:frames2, med_w2, 'LineWidth', 1, 'Color', (color2+white)/2);
            
            subplot(1,3,3); hold on;
            plot(1:frames1, med_diff, 'LineWidth', 1, 'Color', (color3+white)/2);
        end
        
        
        % Per-subject chematic is 2 subplots: 2 pks in one and their difference in the other
        if per_subject_plots
            % Create legends
            trials1 = getSubjectThresholdTrials(subjectId, phases(1));
            leg{1} = [strrep(get_stair_var(phases(1)), '_', ' ') ' (' num2str(sum(trials1)) '/' num2str(length(trials1)) ')'];
            trials2 = getSubjectThresholdTrials(subjectId, phases(1));
            leg{2} = [strrep(get_stair_var(phases(2)), '_', ' ') ' (' num2str(sum(trials2)) '/' num2str(length(trials2)) ')'];
            
            % Plot both PKs
            subplot(1, 2, 1);
            hold on;
            h = boundedline(1:frames1, med_w1', [med_w1-lo_w1; hi_w1-med_w1]', color2, ...
                1:frames2, med_w2', [med_w2-lo_w2; hi_w2-med_w2]', color1, ...
                'alpha');
            errorbar(frames+1, med_b1, med_b1-lo_b1, hi_b1-med_b1, 'LineWidth', 2, 'Color', color2);
            errorbar(frames+1, med_b2, med_b2-lo_b2, hi_b2-med_b2, 'LineWidth', 2, 'Color', color1);
            title([subjectId ' temporal kernels']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 max(frames1, frames2)], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [0 1]);
            xlabel('Time (s)');
            ylabel('Psychophysical Kernel');
            legend(h, leg, 'Location', 'best');
            
            % Plot PK differences
            subplot(1, 2, 2);
            boundedline(1:length(med_diff), med_diff', [med_diff-lo_diff; hi_diff-med_diff]', color3);
            title([subjectId ' kernel difference']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 max(frames1, frames2)], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [-1 0 1]);
            xlabel('Time (s)');
            ylabel('Kernel Difference');
        end
    end
    
    if per_subject_plots
        perSubjectFigs(i).PaperUnits = 'inches';
        perSubjectFigs(i).PaperSize = [8 4];
        perSubjectFigs(i).PaperPosition = [0 0 8 4];
        saveas(perSubjectFigs(i), [subjectId '-PKPlot.png']);
    end
end

if length(subjectIDs) > 1
    figure(combinedFig);
    
    % First, loop over phases and simply concatenate bootstrapped kernels from all subjects together
    for iPhase=length(phases):-1:1
        concatSubjectKernels{iPhase} = vertcat(PerSubjectKernelsByPhase{iPhase, :});
        [~, lo_combo{iPhase}, hi_combo{iPhase}, med_combo{iPhase}] = meanci(concatSubjectKernels{iPhase}, conf);
    end
    
    if length(phases) == 1
        plot(1:frames, med_combo{1}, 'Color', color1, 'LineWidth', 2);
        title([strrep(get_stair_var(phases), '_', ' ') ' combined kernel']);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
    else
        diffAllKernels = concatSubjectKernels{1} - concatSubjectKernels{2};
        [~, lo_diff, hi_diff, med_diff] = meanci(diffAllKernels, conf);
        
        subplot(1, 3, 1);
        hold on;
        plot(1:frames1, med_combo{1}, 'Color', color1, 'LineWidth', 2);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
        title(strrep(get_stair_var(phases(1)), '_', ' '));
        
        subplot(1, 3, 2);
        hold on;
        plot(1:frames2, med_combo{2}, 'Color', color2, 'LineWidth', 2);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames2], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
        title(strrep(get_stair_var(phases(2)), '_', ' '));
        
        subplot(1, 3, 3);
        hold on;
        plot(1:frames1, med_diff, 'Color', color3, 'LineWidth', 2);
        title('Normalized Kernel Difference');
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [-1 0 1]);
        xlabel('Time (s)');
        ylabel('Kernel Differnce');
    end
    
    combinedFig.PaperUnits = 'inches';
    combinedFig.PaperSize = [8 4];
    combinedFig.PaperPosition = [0 0 8 4];
    try
        saveas(combinedFig, [strjoin(subjectIDs, '-') '-PKPlot.png']);
    catch e
        msg = getReport(e);
        warning(msg);
    end
end
end

function stair_var = get_stair_var(phase)
if phase == 0
    stair_var = 'contrast';
elseif phase == 1
    stair_var = 'true_ratio';
elseif phase == 2
    stair_var = 'noise';
else
    error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
end
end

function kernel = normalize(kernel)
kernel = kernel / mean(kernel);
end