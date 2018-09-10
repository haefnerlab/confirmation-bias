function [perSubjectFigs, combinedfig] = DeltaPK(subjectIDs, phases, per_subject_plots, regularize_individuals, method, datadir)
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
if nargin < 4, regularize_individuals = true; end
if nargin < 5, method = 'lr'; end
if nargin < 6, datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

window_low = 0.5;
window_high = 0.7;

    function [median, L, U, trials, frames, variance, true_pk] = getNormalizedSubjectKernel(subjectId, phase, regularize)
        stair_var = get_stair_var(phase);
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        warning off;
        [floor, thresh] = GaborAnalysis.getThresholdWindow(subjectId, phase, window_low, window_high, datadir);
        warning on;
        if phase == 1
            floor = .4;
            thresh = .6;
        end
        trials = SubjectData.(stair_var) <= thresh & SubjectData.(stair_var) >= floor;
        SubjectDataThresh = GaborThresholdTrials(SubjectData, phase, thresh, floor);
        switch lower(method)
            case {'regress', 'logistic', 'lr'}
                if regularize
                    regstring = '-reg';
                    % memo_name = ['PK-xValid-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                    % nFold = 10;
                    % hprs = [0 logspace(-3, 5, 9)];
                    % [hprs, ~] = LoadOrRun(@CustomRegression.xValidatePK, ...
                    %     {SubjectDataThresh.ideal_frame_signals, SubjectDataThresh.choice == +1, hprs, 0, hprs, 1, nFold}, ...
                    %     fullfile(memodir, memo_name));
                    if phase == 1
                        hprs = [0.1 0 10];
                    elseif phase  == 2
                        hprs = [0 0 500];
                    end
                else
                    regstring = '';
                    hprs = [0 0 0];
                end
                memo_name = ['Boot-PK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) regstring '.mat'];
                [~, L, U, median, all_weights] = LoadOrRun(@BootstrapWeightsGabor, ...
                    {SubjectDataThresh, 500, hprs, 0, true}, ...
                    fullfile(memodir, memo_name));
                frames = SubjectData.number_of_images;
                variance = var(all_weights);
            case {'cta'}
                memo_name = ['Boot-CTA-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, L, U, median, all_weights] = LoadOrRun(@BootstrapCTA, ...
                    {SubjectDataThresh, 500}, ...
                    fullfile(memodir, memo_name));
                frames = SubjectData.number_of_images;
                variance = var(all_weights);
            otherwise
                error('Unknown PK method: %s', method);
        end
        if isfield(SubjectData, 'model_pk')
            true_pk = SubjectData.model_pk;
        else
            true_pk = nan;
        end
    end

% CombinedKernelsByPhase{i} contains a the combined kernel mean (weighted by inverse variance pr
% subject) for phase i. Combined kernels are always based on unregularized individual PKs.
PerSubjectKernelsByPhase = cell(length(phases), length(subjectIDs));
CombinedKernelsByPhase = cell(size(phases));
CombinedNormalizerByPhase = cell(size(phases));
CombinedIdentifier = strjoin(subjectIDs, '-');

perSubjectFigs = [];
for i=1:length(subjectIDs)
    if per_subject_plots
        perSubjectFigs(i) = figure;
        hold on;
    end
    subjectId = subjectIDs{i};
    if length(phases) == 1
        [M, L, U, trials, frames, variance, true_pk] = getNormalizedSubjectKernel(subjectId, phases, false);
        if per_subject_plots
            boundedline(1:frames, M(1:frames)', [U(1:frames)-M(1:frames); M(1:frames)-L(1:frames)]');
            errorbar(frames+1, M(end), M(end)-L(end), U(end)-M(end), 'LineWidth', 2, 'Color', 'r');
            title([strrep(get_stair_var(phases), '_', ' ') 'temporal kernel (' num2str(sum(trials)) '/' num2str(length(trials)) ')']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 frames], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [0 1]);
            xlabel('Time (s)');
            ylabel('Psychophysical Kernel');
            axis tight;
            if ~isnan(true_pk)
                plot(1:frames, true_pk, 'Color', 'r', 'LineWidth', 2);
            end
        end
        if i == 1
            CombinedKernelsByPhase{1} = M ./ variance;
            CombinedNormalizerByPhase{1} = 1 ./ variance;
        else
            CombinedKernelsByPhase{1} = CombinedKernelsByPhase{1} + M ./ variance;
            CombinedNormalizerByPhase{1} = CombinedNormalizerByPhase{1} + 1 ./ variance;
        end
        
        if regularize_individuals
            [PerSubjectKernelsByPhase{1, i}, ~, ~, ~, ~, ~, ~] = ...
                getNormalizedSubjectKernel(subjectId, phases, true);
        else
            PerSubjectKernelsByPhase{1, i} = M;
        end
    else
        % 2 subplots: 2 pks in one and their difference in the other
        [M1, L1, U1, trials1, frames1, v1, true_pk1] = getNormalizedSubjectKernel(subjectId, phases(1), false);
        leg{1} = [strrep(get_stair_var(phases(1)), '_', ' ') ' (' num2str(sum(trials1)) '/' num2str(length(trials1)) ')'];
        
        [M2, L2, U2, trials2, frames2, v2, true_pk2] = getNormalizedSubjectKernel(subjectId, phases(2), false);
        leg{2} = [strrep(get_stair_var(phases(2)), '_', ' ') ' (' num2str(sum(trials2)) '/' num2str(length(trials2)) ')'];
        
        if per_subject_plots
            subplot(1, 2, 1);
            hold on;
            h = boundedline(1:frames1, M1(1:frames1)', [U1(1:frames1)-M1(1:frames1); M1(1:frames1)-L1(1:frames1)]', 'b', ...
                1:frames2, M2(1:frames2)', [U2(1:frames2)-M2(1:frames2); M2(1:frames2)-L2(1:frames2)]', 'r', ...
                'alpha');
            title([subjectId ' temporal kernels']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [0 1]);
            xlabel('Time (s)');
            ylabel('Psychophysical Kernel');
            legend(h, leg, 'Location', 'best');
            if ~isnan(true_pk1)
                plot(1:frames1, true_pk1, 'Color', 'b', 'LineWidth', 2);
                plot(1:frames2, true_pk2, 'Color', 'r', 'LineWidth', 2);
            end
            
            subplot(1, 2, 2);
            assert(frames1 == frames2, 'Cannot subtract PKs with different # frames');
            Md = M1 - M2;
            Ld = L1 - U2;
            Ud = U1 - L2;
            boundedline(1:frames1, Md(1:frames1)', [Ud(1:frames1)-Md(1:frames1); Md(1:frames1)-Ld(1:frames1)]', 'k');
            title([subjectId ' kernel difference']);
            xlim([-inf, inf]);
            set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
            set(gca, 'YTick', [-1 0 1]);
            xlabel('Time (s)');
            ylabel('Kernel Difference');
        end
        
        if i == 1
            CombinedKernelsByPhase{1} = M1 ./ v1;
            CombinedNormalizerByPhase{1} = 1 ./ v1;
            CombinedKernelsByPhase{2} = M2 ./ v2;
            CombinedNormalizerByPhase{2} = 1 ./ v2;
        else
            CombinedKernelsByPhase{1} = CombinedKernelsByPhase{1} + M1 ./ v1;
            CombinedNormalizerByPhase{1} = CombinedNormalizerByPhase{1} + 1 ./ v1;
            CombinedKernelsByPhase{2} = CombinedKernelsByPhase{2} + M2 ./ v2;
            CombinedNormalizerByPhase{2} = CombinedNormalizerByPhase{2} + 1 ./ v2;
        end
        
        
        if regularize_individuals
            [PerSubjectKernelsByPhase{1, i}, ~, ~, ~, ~, ~, ~] = ...
                getNormalizedSubjectKernel(subjectId, phases(1), true);
            [PerSubjectKernelsByPhase{2, i}, ~, ~, ~, ~, ~, ~] = ...
                getNormalizedSubjectKernel(subjectId, phases(2), true);
        else
            PerSubjectKernelsByPhase{1, i} = M1;
            PerSubjectKernelsByPhase{2, i} = M2;
        end
    end
    
    if per_subject_plots
        perSubjectFigs(i).PaperUnits = 'inches';
        perSubjectFigs(i).PaperSize = [8 4];
        perSubjectFigs(i).PaperPosition = [0 0 8 4];
        saveas(perSubjectFigs(i), [subjectId '-PKPlot.png']);
    end
end

combinedfig = -1;
if length(subjectIDs) > 1
    combinedfig = figure;
    hold on;
    
    if length(phases) == 1
        % Normalization
        CombinedKernelsByPhase{1} = CombinedKernelsByPhase{1} ./ CombinedNormalizerByPhase{1};
        % Plot
        plot(1:frames, CombinedKernelsByPhase{1}(1:frames));
        title([strrep(get_stair_var(phases), '_', ' ') ' combined kernel']);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
    else
        colors = [0 0 1; 1 0 0; 0 0 0];
        % Normalization
        CombinedKernelsByPhase{1} = CombinedKernelsByPhase{1} ./ CombinedNormalizerByPhase{1};
        CombinedKernelsByPhase{2} = CombinedKernelsByPhase{2} ./ CombinedNormalizerByPhase{2};
        % 3 subplots: each condition per-subject and mean PKs and their difference
        subplot(1, 3, 1);
        hold on;
        for i=1:length(subjectIDs)
            plot(1:frames1, PerSubjectKernelsByPhase{1, i}(1:frames1), 'Color', colors(1, :));
        end
        plot(1:frames1, CombinedKernelsByPhase{1}(1:frames1), 'Color', colors(1, :), 'LineWidth', 2);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
        title(get_stair_var(phases(1)));
        
        subplot(1, 3, 2);
        hold on;
        for i=1:length(subjectIDs)
            plot(1:frames2, PerSubjectKernelsByPhase{2, i}(1:frames2), 'Color', colors(2, :));
        end
        plot(1:frames2, CombinedKernelsByPhase{2}(1:frames2), 'Color', colors(2, :), 'LineWidth', 2);
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [0 1]);
        xlabel('Time (s)');
        ylabel('Psychophysical Kernel');
        title(get_stair_var(phases(2)));
        
        subplot(1, 3, 3);
        hold on;
        assert(frames1 == frames2, 'Cannot subtract PKs with different # frames');
        for i=1:length(subjectIDs)
            diff = normalize(PerSubjectKernelsByPhase{1, i}(1:frames1)) - ...
                normalize(PerSubjectKernelsByPhase{2, i}(1:frames2));
            plot(1:frames1, diff, 'Color', colors(3, :));
        end
        combo_diff = normalize(CombinedKernelsByPhase{1}(1:frames1)) - ...
            normalize(CombinedKernelsByPhase{2}(1:frames2));
        plot(1:frames1, combo_diff, 'Color', colors(3, :), 'LineWidth', 2);
        title('Normalized Kernel Difference');
        xlim([-inf, inf]);
        set(gca, 'XTick', [1 frames1], 'XTickLabel', [0 1]);
        set(gca, 'YTick', [-1 0 1]);
        xlabel('Time (s)');
        ylabel('Kernel Differnce');
    end
    
    combinedfig.PaperUnits = 'inches';
    combinedfig.PaperSize = [8 4];
    combinedfig.PaperPosition = [0 0 8 4];
    try
        saveas(combinedfig, [CombinedIdentifier '-PKPlot.png']);
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