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
        ax = subplot(nS, nP, sub2ind([nS nP], i, j));
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
            case 'template'
            case 'pk'
        end
    end
end

% combined = figure();
end