function [grid, combined] = PlotGrid(subjectIDs, thresholds, phase, plot_types, datadir)
%GABORANALYSIS.PLOTGRID creates two figures: one a grid of subplots with
%each subject as a row and each plot type as a column. The second shows the
%combined PK for all subjects. Returns figure handles. Example:
%
% GABORANALYSIS.PLOTGRID(subjectIDs, thresholds, phase, plot_types, datadir)
%
% Inputs:
% - subjectIDs: a cell array of subject IDs (strings)
% - thresholds: an array the same size as subjectIDs with each subject's
%   threshold, or a 2d array (subjects x 2) with [lo hi] ranges to keep.
% - phase: 0 for contrast, 1 for ratio
% - plot_types: (optional) cell array of plots to show. options are 'staircase', 'rt', 'pm', 'pk', 'sd'
% - datadir: (optional) override the default place to look for data files.

if nargin < 4, plot_types = {'staircase', 'pm', 'pk'}; end
if nargin < 5, datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

if phase == 0
    stair_var = 'contrast';
elseif phase == 1
    stair_var = 'true_ratio';
elseif phase == 2
    stair_var = 'noise';
else
    error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
end

if size(thresholds, 2) == 1
    thresholds = [zeros(length(subjectIDs), 1), thresholds(:)];
end

nS = length(subjectIDs);
nP = length(plot_types);

window_low = 0.5;
window_high = 0.75;

    function [floor, thresh] = getThresholdWrapper(subjectId)
        s_idx = strcmpi(subjectId, subjectIDs);
        if isinf(thresholds(s_idx, 2))
            [floor, thresh] = GaborAnalysis.getThresholdWindow(subjectId, phase, window_low, window_high, datadir);
        else
            thresh = thresholds(s_idx, 2);
            if phase == 1
                floor = 1 - thresh;
            else
                floor = thresholds(s_idx, 1);
            end
        end
    end

grid = figure();
for i=1:nS
    s = subjectIDs{i};
    SubjectData = LoadOrRun(@LoadAllSubjectData, ...
        {s, phase, datadir}, fullfile(catdir, [s '-' stair_var '.mat']));
    [floor, thresh] = getThresholdWrapper(s);
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
                [fit_result, uniq_vals, yvals, stderrs] = LoadOrRun(@GaborPsychometric, ...
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
                    plotOptions.xLabel = 'True Ratio';
                    plotOptions.yLabel = 'Percent Chose Left';
                    
                    % Plot PM curve and data.
                    plotPsych(fit_result, plotOptions);
                    
                    % Plot data.
                    errorbar(uniq_vals, yvals, stderrs, 'bs');
                    ys = get(gca, 'YLim');
                    plot([floor floor], ys, '--r');
                    plot([thresh thresh], ys, '--r');
                elseif phase == 2
                    % Add remaining plot options.
                    plotOptions.xLabel = 'Noise (\kappa)';
                    plotOptions.yLabel = 'Percent Correct';
                    
                    % Plot PM curve and data.
                    plotPsych(fit_result, plotOptions);
                    
                    % Plot data.
                    errorbar(uniq_vals, yvals, stderrs, 'bs');
                    ys = get(gca, 'YLim');
                    plot([floor floor], ys, '--r');
                    plot([thresh thresh], ys, '--r');
                end
                title('Psychometric curve');
            case {'pk', 'pk-lr'}
                SubjectDataThresh = GaborThresholdTrials(...
                    SubjectData, phase, thresh, floor);
                memo_name = ['Boot-PK-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [M, L, U, ~] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, 500}, ...
                    fullfile(memodir, memo_name));
                frames = SubjectData.number_of_images;
                boundedline(1:frames, M(1:frames)', [U(1:frames)-M(1:frames); M(1:frames)-L(1:frames)]');
                errorbar(frames+1, M(end), M(end)-L(end), U(end)-M(end), 'LineWidth', 2, 'Color', 'r');
                title(['temporal kernel (LR)']);
                set(gca, 'XAxisLocation', 'origin');
                set(gca, 'XTick', [1 frames]);
                axis tight;
            case {'cta', 'pk-cta'}
                SubjectDataThresh = GaborThresholdTrials(...
                    SubjectData, phase, thresh, floor);
                memo_name = ['Boot-CTA-' stair_var '-' s '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [M, L, U, ~] = LoadOrRun(@BootstrapCTA, ...
                    {SubjectDataThresh, 500}, ...
                    fullfile(memodir, memo_name));
                frames = SubjectData.number_of_images;
                boundedline(1:frames, M(1:frames)', [U(1:frames)-M(1:frames); M(1:frames)-L(1:frames)]');
                errorbar(frames+1, M(end), M(end)-L(end), U(end)-M(end), 'LineWidth', 2, 'Color', 'r');
                title(['temporal kernel (CTA)']);
                set(gca, 'XAxisLocation', 'origin');
                set(gca, 'XTick', [1 frames]);
                axis tight;
            case 'streaks'
                correct_streaks = [];
                incorrect_streaks = [];
                mode = 0;
                current_streak = 0;
                for t=1:SubjectData.current_trial
                    if SubjectData.accuracy(t)
                        if mode == +1
                            current_streak = current_streak + 1;
                        elseif mode == -1
                            incorrect_streaks(end+1) = current_streak;
                            current_streak = 1;
                        end
                        mode = +1;
                    else
                        if mode == -1
                            current_streak = current_streak + 1;
                        elseif mode == +1
                            correct_streaks(end+1) = current_streak;
                            current_streak = 1;
                        end
                        mode = -1;
                    end
                end
                mu = mean(correct_streaks);
                xs = linspace(1, max(correct_streaks));
                correct_exp_pdf = exppdf(xs, mu);
                [ux, uc] = unique_counts(correct_streaks);
                bar(ux, uc/sum(uc), .5, 'k');
                
                mu = mean(incorrect_streaks);
                incorrect_exp_pdf = exppdf(xs, mu);
                [ux, uc] = unique_counts(incorrect_streaks);
                bar(ux+.5, uc/sum(uc), .5, 'r');
                
                plot(xs, correct_exp_pdf, '--k');
                plot(xs, incorrect_exp_pdf, '--r');
                
                legend('correct', 'incorrect');
                xlabel('# consecutive');
                title('histograms of streakiness');
            case 'ideal-bias'
                unq_signals = unique(SubjectData.(stair_var));
                bias = zeros(3, length(unq_signals));
                for ui=1:length(unq_signals)
                    trials_at = SubjectData.(stair_var) == unq_signals(ui);
                    ntrials = sum(trials_at);
                    left = sum(SubjectData.correct_answer(trials_at));
                    [bias(1, ui), bias(2:3, ui)] = binofit(left, ntrials);
                end
                bar(1:length(unq_signals), bias(1, :));
                errorbar(1:length(unq_signals), bias(1, :), bias(1, :) - bias(2, :), bias(3, :) - bias(1, :), '.');
            case 'ideal-signals'
                unq_signals = unique(SubjectData.(stair_var));
                for ui=1:length(unq_signals)
                    trials_at = SubjectData.(stair_var) == unq_signals(ui);
                    frame_signals = SubjectData.ideal_frame_signals(trials_at, :);
                    ksdensity(frame_signals(:));
                end
                xlabel('Per Frame Signal');
                legend(arrayfun(@(s) [stair_var '=' num2str(s)], unq_signals, 'UniformOutput', false));
                axis tight;
        end
        drawnow;
    end
end
end