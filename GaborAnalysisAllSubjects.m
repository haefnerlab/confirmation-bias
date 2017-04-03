function [grid, combined] = GaborAnalysisAllSubjects(subjectIDs, thresholds, phase, plot_types, ideal_template, datadir)
%GABORANALYSISALLSUBJECTS creates two figures: one a grid of subplots with
%each subject as a row and each plot type as a column. The second shows the
%combined PK for all subjects. Returns figure handles. Example:
%
% GABORANALYSISALLSUBJECTS(subjectIDs, thresholds, phase, plot_types, ideal_template, datadir)
%
% Inputs:
% - subjectIDs: a cell array of subject IDs (strings)
% - thresholds: an array the same size as subjectIDs with each subject's
%   threshold, or a 2d array (subjects x 2) with [lo hi] ranges to keep.
% - phase: 0 for contrast, 1 for ratio
% - plot_types: (optional) cell array of plots to show. options are 'staircase', 'rt', 'pm', 'template', 'pk', 'sd'
% - ideal_template: (optional) if true, uses ideal observer spatial template.
% - datadir: (optional) override the default place to look for data files.

if nargin < 4, plot_types = {'staircase', 'pm', 'template', 'pk'}; end
if nargin < 5, ideal_template = true; end
if nargin < 6, datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

if phase == 0
    stair_var = 'contrast';
elseif phase == 1
    stair_var = 'true_ratio';
else
    error('Expected phase 0 for Contrast or 1 for Ratio');
end

if numel(thresholds) == length(subjectIDs)
    thresholds = [zeros(length(subjectIDs), 1), thresholds(:)];
end

nS = length(subjectIDs);
nP = length(plot_types);

    function [floor, thresh] = getThresholdWindow(subjectID, perf_lo, perf_hi)
        subjectIdx = strcmpi(subjectID, subjectIDs);
        floor = thresholds(subjectIdx, 1);
        thresh = thresholds(subjectIdx, 2);
        if isinf(thresh)
            [LocalSubjectData, ~] = ...
                LoadOrRun(@LoadAllSubjectData, {subjectID, phase, datadir}, fullfile(catdir, [subjectID '-' stair_var '.mat']));
            % Use PM fit to get floor and threshold
            local_fit_result = LoadOrRun(@GaborPsychometric, ...
                {LocalSubjectData, phase}, ...
                fullfile(memodir, ['PM-' stair_var '-' subjectID '.mat']));
            floor = getThreshold(local_fit_result, perf_lo, false);
            thresh = getThreshold(local_fit_result, perf_hi, false);
        
            % Adjust from '#clicks' to threshold
            if phase == 1
                thresh = thresh / 10;
                floor = 1 - thresh;
            end
        end
    end

grid = figure();
for i=1:nS
    s = subjectIDs{i};
    [SubjectData, stim_images] = ...
        LoadOrRun(@LoadAllSubjectData, {s, phase, datadir}, fullfile(catdir, [s '-' stair_var '.mat']));
    [floor, thresh] = getThresholdWindow(s, .6, .75);
    trials = SubjectData.(stair_var) <= thresh & SubjectData.(stair_var) >= floor;
    
    for j=1:length(plot_types)
        ax = subplot(nS, nP, (i-1)*nP + j);
        hold on;
        switch lower(plot_types{j})
            case 'staircase'
                plot(SubjectData.(stair_var), '-k');
                plot(find(trials), SubjectData.(stair_var)(trials), 'ob');
                line([1 length(trials)], [thresh thresh], 'LineStyle', '--', 'Color', 'r');
                if floor > 0
                    line([1 length(trials)], [floor floor], 'LineStyle', '--', 'Color', 'r');
                end
                ylabel(stair_var);
                title([num2str(sum(trials)) '/' num2str(length(trials)) ' trials']);
            case 'rt'
                plot(SubjectData.reaction_time);
                title('reaction time');
            case 'sd'
                [pCL, pCR, pWL, pWR] = Serial_Dependencies(SubjectData);
                bar([pCL, pCR, pWL, pWR]);
                xlabel('Previous Trial');
                ylabel('Prob(choose left)');
                set(ax, 'XTick', 1:4);
                set(ax, 'XTickLabel', {'Left+Correct', 'Right+Correct', 'Left+Wrong', 'Right+Wrong'});
                if exist('xtickangle', 'file'), xtickangle(ax, 25); end
                title('serial dependencies');
            case 'pm'
                % Get PM fit.
                fit_result = LoadOrRun(@GaborPsychometric, ...
                    {SubjectData, phase}, ...
                    fullfile(memodir, ['PM-' stair_var '-' s '.mat']));
                
                % Construct options for psignifit plotting.
                plotOptions = struct;
                plotOptions.plotData       = false;
                plotOptions.plotAsymptote  = false;
                plotOptions.plotThresh     = false;
                plotOptions.CIthresh       = false;
                
                if phase == 0
                    % Add remaining plot options.
                    plotOptions.xLabel = 'Contrast Level';
                    plotOptions.yLabel = 'Percent Correct';
                    
                    % Plot PM curve and data.
                    plotPsych(fit_result, plotOptions);
                    
                    % Bin the data further for visualization
                    log_contrast = log(SubjectData.contrast);
                    bin_edges = linspace(min(log_contrast), max(log_contrast), 11);
                    bin_halfwidth = (bin_edges(2) - bin_edges(1)) / 2;
                    bin_centers = bin_edges(1:end-1) + bin_halfwidth;
                    means = zeros(size(bin_centers));
                    stderrs = zeros(size(bin_centers));
                    for b=1:length(bin_centers)
                        % Select all points for which bin i is closest.
                        bin_dists = abs(log_contrast - bin_centers(b));
                        indices = bin_dists <= bin_halfwidth;
                        means(b) = mean(SubjectData.accuracy(indices));
                        stderrs(b) = std(SubjectData.accuracy(indices)) / sqrt(sum(indices));
                    end
                    errorbar(exp(bin_centers), means, stderrs, 'bs');
                    ys = get(gca, 'YLim');
                    plot([floor floor], ys, '--r');
                    plot([thresh thresh], ys, '--r');
                elseif phase == 1
                    % Add remaining plot options.
                    plotOptions.xLabel = 'True # Left Frames';
                    plotOptions.yLabel = 'Percent Chose Left';
                    
                    % Plot PM curve and data.
                    plotPsych(fit_result, plotOptions);
                    
                    % Plot data.
                    true_ratio = sum(SubjectData.order_of_orientations, 2);
                    uniq_vals = unique(true_ratio);
                    
                    % The next line exploits the fact that 'Left' is coded as 1 and 'right'
                    % as 0, so mean() returns the fraction of left choices.
                    yvals = arrayfun(@(u) mean(SubjectData.choice(true_ratio == u)), uniq_vals);
                    num_trials_at_vals = arrayfun(@(u) sum(true_ratio == u), uniq_vals);
                    stderrs = arrayfun(@(u) std(SubjectData.choice(true_ratio == u)), uniq_vals) ./ sqrt(num_trials_at_vals);
                    errorbar(uniq_vals, yvals, stderrs, 'bs');
                    ys = get(gca, 'YLim');
                    plot(10*[floor floor], ys, '--r');
                    plot(10*[thresh thresh], ys, '--r');
                end
                title('Psychometric curve');
                
                % Display floor+threshold values
                [floor, thresh] = getThresholdWindow(s, .6, .75);
            case 'template'
                [SubjectDataThresh, images_thresh] = GaborThresholdTrials(...
                    SubjectData, stim_images, phase, thresh, floor);
                if ideal_template
                    memo_name = ['Boot-PK-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                else
                    memo_name = ['Boot-SpPK-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                end
                [M, ~, ~] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, images_thresh, 500, ideal_template}, ...
                    fullfile(memodir, memo_name));
                [~, ~, h, w] = size(stim_images);
                template = reshape(M(1:h*w), [h w]);
                imagesc(template);
                axis image;
                colorbar;
                title('spatial kernel');
            case 'pk'
                [SubjectDataThresh, images_thresh] = GaborThresholdTrials(...
                    SubjectData, stim_images, phase, thresh, floor);
                if ideal_template
                    memo_name = ['Boot-PK-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                else
                    memo_name = ['Boot-SpPK-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                end
                [M, L, U] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, images_thresh, 500, ideal_template}, ...
                    fullfile(memodir, memo_name));
                [~, frames, h, w] = size(stim_images);
                time_idxs = w*h+1:length(M)-1;
                boundedline(1:frames, M(time_idxs)', [U(time_idxs)-M(time_idxs); M(time_idxs)-L(time_idxs)]');
                title('temporal kernel');
        end
    end
end

if length(subjectIDs) > 1 && any(strcmpi('pk', plot_types))
    combined = figure();
    all_data = [];
    all_imgs = [];
    all_ids = '';
    for i=1:nS
        s = subjectIDs{i};
        [SubjectData, stim_images] = ...
            LoadOrRun(@LoadAllSubjectData, {s, phase, datadir}, fullfile(catdir, [s '-' stair_var '.mat']));
        [floor, thresh] = getThresholdWindow(s, .6, .75);
        all_ids = strcat(all_ids, ['-' s '-' num2str(thresh) '-' num2str(floor)]);
        [this_subject, this_imgs] = GaborThresholdTrials(SubjectData, stim_images, phase, thresh, floor);
        if isempty(all_data)
            all_data = this_subject;
            all_imgs = this_imgs;
        else
            [all_data, all_imgs] = ConcatGaborData(all_data, all_imgs, this_subject, this_imgs);
        end
    end
    [M, L, U] = LoadOrRun(@BootstrapWeightsGabor, ...
        {all_data, all_imgs, 500, ideal_template}, ...
        fullfile(memodir, ['combo-PK-' stair_var all_ids '.mat']));
    [~, frames, h, w] = size(all_imgs);
    time_idxs = w*h+1:length(M)-1;
    boundedline(1:frames, M(time_idxs)', [U(time_idxs)-M(time_idxs); M(time_idxs)-L(time_idxs)]');
    title('combined temporal kernel');
    
elseif nargout > 1
    error('Need ''pk'' plot type and at least 2 subjects to make a combined plot');
end
end