function [slopes, pvalues] = DeltaSlopeStatistics(subjectIds, phases, type, nboot, datadir)
%DELTASLOPESTATISTICS use linear PK fit to each of the two phases, bootstrapped, and return (1) the
%[num subjects x 2 x bootstraps] array of PK slopes and (2) the [num subjects] array of p-values for
%there being a significant difference between the two phases.

if ~exist('nboot', 'var'), nboot = 10000; end
if ~exist('datadir', 'var'), datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

slopes = nan(length(subjectIds), 2, nboot);
pvalues = nan(length(subjectIds), 1);

window_low = 0.5;
window_high = 0.7;

    function [floor, thresh] = getThresholdWrapper(subjectId, phase)
        if phase == 1
            floor = .4;
            thresh = .6;
        else
            [floor, thresh] = GaborAnalysis.getThresholdWindow(subjectId, phase, window_low, window_high, datadir);
        end
    end

    function [stair_var, stair_var_name] = getStairVar(phase)
        if phase == 0
            stair_var = 'contrast';
            stair_var_name = 'contrast';
        elseif phase == 1
            stair_var = 'true_ratio';
            stair_var_name = 'ratio';
        elseif phase == 2
            stair_var = 'noise';
            stair_var_name = 'noise';
        else
            error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
        end
    end

    function slopes = getSlopesForPhase(subjectId, phase)
        stair_var = getStairVar(phase);
        
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        [floor, thresh] = getThresholdWrapper(subjectId, phase);
        SubjectDataThresh = GaborThresholdTrials(SubjectData, phase, thresh, floor);
        
        switch type
            case 'linear'
                memo_name = ['Boot-LinPK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, ~, fits] = LoadOrRun(@BootstrapLinearPKFit, {SubjectDataThresh, nboot, 0, true}, ...
                    fullfile(memodir, memo_name));
                slopes = fits(:, 1);
            case 'exponential'
                memo_name = ['Boot-ExpPK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
                [~, ~, ~, ~, fits] = LoadOrRun(@BootstrapExponentialWeightsGabor, {SubjectDataThresh, nboot, 0, true}, ...
                    fullfile(memodir, memo_name));
                slopes = fits(:, 2);
        end
    end

    function name = cleanName(subjectId)
        st = regexpi(subjectId, 'subject\d+');
        name = subjectId(st:end);
    end

switch type
    case 'linear'
        xlab = 'linear slope';
        range = [-.2 .2];
    case 'exponential'
        xlab = '\beta';
        range = [-1 1];
end

figure;
edges = range(1):.005:range(2);
edgediff = edges(2) - edges(1);
all_slopes_diff = zeros(length(subjectIds), nboot);
all_slopes_diff_medians = zeros(length(subjectIds), 1);
for s_idx=1:length(subjectIds)
    %% Compute fits for subject
    slopes1 = getSlopesForPhase(subjectIds{s_idx}, phases(1));
    slopes(s_idx, 1, :) = slopes1;
    
    slopes2 = getSlopesForPhase(subjectIds{s_idx}, phases(2));
    slopes(s_idx, 2, :) = slopes2;
    
    %% Take difference of slopes and compute pvalue
    slopes_diff = slopes1 - slopes2;
    pvalues(s_idx) = (1 + sum(slopes_diff < 0)) / (1 + nboot);
    
    all_slopes_diff(s_idx, :) = slopes_diff;
    all_slopes_diff_medians(s_idx) = median(slopes_diff);
    
    %% Plot histograms
    subplotsquare(length(subjectIds), s_idx);
    cla; hold on;
    if pvalues(s_idx) == 1 / (nboot + 1)
        p_string = ['p < ' num2str(1/nboot)];
    elseif pvalues(s_idx) < .01
        p_string = 'p < .01';
    else
        p_string = ['p = ' num2str(pvalues(s_idx))];
    end
    title(sprintf('%s (%s)', cleanName(subjectIds{s_idx}), p_string));
    
    hist1 = histcounts(slopes1, edges);
    hist2 = histcounts(slopes2, edges);
    bar(edges(1:end-1), hist1, .5, 'b', 'EdgeColor', 'none');
    bar(edges(1:end-1)+.5*edgediff, hist2, .5, 'r', 'EdgeColor', 'none');
    xlabel(['PK slope (' xlab ')']);
    ylabel('bootstrap histogram');
    set(gca, 'YTick', 0);
    xlim(range);
    
    if s_idx == 1
        [~, label1] = getStairVar(phases(1));
        [~, label2] = getStairVar(phases(2));
        legend(label1, label2);
    end
end

%% Histogram ALL the subjects
figure;
hold on;
binedges = linspace(range(1), range(2), 120);
binwidth = (binedges(2) - binedges(1));
bincenters = binedges(1:end-1) + binwidth / 2;
counts = histcounts(all_slopes_diff(:), binedges);
% --- Each subject as a different color ---
% per_counts = zeros(length(subjectIds), length(bincenters));
% for s_idx=1:length(subjectIds)
%     per_counts(s_idx, :) = histcounts(all_slopes_diff(s_idx, :), binedges);
% end
% bar(bincenters, per_counts', 'stacked');
% % --- All subjects plus ksdensity over medians ---
% bar(bincenters, counts / sum(counts) / binwidth); % All subjects x all bootstraps combined
% [density, dens_pts] = ksdensity(all_slopes_diff_medians);
% plot(dens_pts, density, '-r', 'LineWidth', 2);
% --- KS density of each subject as a different color ---
for s_idx=1:length(subjectIds)
    [density, dens_pts] = ksdensity(all_slopes_diff(s_idx, :));
    plot(dens_pts, density, 'LineWidth', 2);
end
xlabel(['\Delta PK Slope (' xlab ')']);
xlim(range);

%% Violin plots

% Take mean over bootstrap samples
subj_mean_slopes = mean(slopes, 3);

figure;
violins = violinplot(subj_mean_slopes, {label1, label2}, 'ShowData', false, 'ShowMean', false, 'ShowNotches', false, 'ViolinAlpha', 1);
violins(1).ViolinColor = SamplingModel.betacolor(1, -1, 1);
delete(violins(1).BoxPlot);
delete(violins(1).MedianPlot);
violins(2).ViolinColor = SamplingModel.betacolor(-1, -1, 1);
delete(violins(2).BoxPlot);
delete(violins(2).MedianPlot);

hold on;
for iSubject=1:length(subjectIds)
    plot([1 2], subj_mean_slopes(iSubject, :), '-ok', 'MarkerFaceColor', 'k');
end

end