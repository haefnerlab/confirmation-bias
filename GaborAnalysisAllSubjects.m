function [grid, combined] = GaborAnalysisAllSubjects(subjectIDs, thresholds, phase, plot_types, ideal_template, datadir)
%GABORANALYSISALLSUBJECTS creates two figures: one a grid of subplots with
%each subject as a row and each plot type as a column. The second shows the
%combined PK for all subjects. Returns figure handles. Example:
%
% GABORANALYSISALLSUBJECTS(subjectIDs, thresholds, phase, plot_types, ideal_template, datadir)
%
% Inputs:
% - subjectIDs: a cell array of subject IDs (strings)
% - thresholds: an array the same size as subjectIDs with each subject's threshold.
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
    stair_var = 'ratio';
else
    error('Expected phase 0 for Contrast or 1 for Ratio');
end

nS = length(subjectIDs);
nP = length(plot_types);

grid = figure();
for i=1:nS
    s = subjectIDs{i};
    [SubjectData, stim_images] = ...
        LoadOrRun(@LoadAllSubjectData, {s, phase, datadir}, fullfile(catdir, [s '.mat']));
    thresh = thresholds(i);
    trials = SubjectData.(stair_var) <= thresh;
    
    for j=1:length(plot_types)
        ax = subplot(nS, nP, (i-1)*nP + j);
        hold on;
        switch lower(plot_types{j})
            case 'staircase'
                plot(SubjectData.(stair_var), '-k');
                plot(find(trials), SubjectData.(stair_var)(trials), 'ob');
                line([1 length(trials)], [thresh thresh], 'LineStyle', '--', 'Color', 'r');
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
                xtickangle(ax, 25);
                title('serial dependencies');
            case 'pm'
                % Construct options for psignifit function and subsequent plotting.
                options = struct;
                plotOptions = struct;
                plotOptions.plotData       = false;
                plotOptions.plotAsymptote  = false;
                plotOptions.plotThresh     = false;
                plotOptions.CIthresh       = false;
                
                if phase == 0
                    options.sigmoidName  = 'weibull';
                    options.expType      = '2AFC';
                    % Count how often subject was correct at each contrast
                    % value.
                    uniq_vals = unique(SubjectData.contrast);
                    yvals = arrayfun(@(u) mean(SubjectData.accuracy(SubjectData.contrast == u)), uniq_vals);
                    num_trials_at_vals = arrayfun(@(u) sum(SubjectData.contrast == u), uniq_vals);
                    
                    % Add remaining plot options.
                    plotOptions.xLabel = 'Contrast Level';
                    plotOptions.yLabel = 'Percent Correct';
                    options.estimateType = 'MLE';
                    options.nblocks      = length(uniq_vals);
                    options.threshPC     = .7;
                    options.betaPrior    = 10;
                    
                    % Run PM fitting.
                    result = LoadOrRun(@psignifit, ...
                        {[uniq_vals(:) yvals(:) num_trials_at_vals(:)], options}, ...
                        fullfile(memodir, ['PM-' stair_var '-' s '.mat']));
                    
                    % Plot PM curve and data.
                    plotPsych(result, plotOptions);
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
                    errorbar(exp(bin_centers), means, stderrs, 'bs');  % Black line
                elseif phase == 1
                    options.sigmoidName  = 'norm';
                    options.expType      = 'YesNo';
                    % Count how often subject chose left at each ratio
                    % value.
                    true_ratio = sum(SubjectData.order_of_orientations, 2);
                    uniq_vals = unique(true_ratio);
                    % The next line exploits the fact that 'Left' is coded
                    % as 1 and 'right' as 0, so mean() returns the fraction
                    % of left choices.
                    yvals = arrayfun(@(u) mean(SubjectData.choice(true_ratio == u)), uniq_vals);
                    num_trials_at_vals = arrayfun(@(u) sum(true_ratio == u), uniq_vals);
                    stderr = arrayfun(@(u) std(SubjectData.choice(true_ratio == u)), uniq_vals) ./ sqrt(num_trials_at_vals);
                    
                    % Add remaining plot options.
                    plotOptions.xLabel = 'True # Left Frames';
                    plotOptions.yLabel = 'Percent Chose Left';
                    options.estimateType = 'MLE';
                    options.nblocks      = length(uniq_vals);
                    options.threshPC     = .7;
                    options.betaPrior    = 10;
                    
                    % Run PM fitting.
                    result = LoadOrRun(@psignifit, ...
                        {[uniq_vals(:) yvals(:) num_trials_at_vals(:)], options}, ...
                        fullfile(memodir, ['PM-' stair_var '-' s '.mat']));
                    
                    % Plot PM curve and data.
                    plotPsych(result, plotOptions);
                    errorbar(uniq_vals, yvals, stderr, 'bs');  % Black line
                end
                title('Psychometric curve');
            case 'template'
            case 'pk'
        end
    end
end

% combined = figure();
end