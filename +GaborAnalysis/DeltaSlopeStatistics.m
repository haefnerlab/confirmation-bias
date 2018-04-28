function [slopes, pvalues] = DeltaSlopeStatistics(subjectIds, phases, nboot, datadir)
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
window_high = 0.75;

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

    function fits = doFitsForPhase(subjectId, phase)
    stair_var = getStairVar(phase);
    
    SubjectData = LoadOrRun(@LoadAllSubjectData, ...
        {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
    [floor, thresh] = getThresholdWrapper(subjectId, phase);
    SubjectDataThresh = GaborThresholdTrials(SubjectData, phase, thresh, floor);
    
    memo_name = ['Boot-LinPK-ideal-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [~, ~, ~, ~, ~, fits] = LoadOrRun(@BootstrapLinearPKFit, {SubjectDataThresh, nboot}, ...
        fullfile(memodir, memo_name));
    end

    function name = cleanName(subjectId)
    st = regexpi(subjectId, 'subject\d+');
    name = subjectId(st:end);
    end

edges = -.2:.005:.2;
edgediff = edges(2) - edges(1);
for s_idx=1:length(subjectIds)
    %% Compute fits for subject
    fits1 = doFitsForPhase(subjectIds{s_idx}, phases(1));
    slopes1 = fits1(:, 1);
    slopes(s_idx, 1, :) = slopes1;
    
    fits2 = doFitsForPhase(subjectIds{s_idx}, phases(2));
    slopes2 = fits2(:, 1);
    slopes(s_idx, 2, :) = slopes2;
    
    %% Take difference of slopes and compute pvalue
    slopes_diff = slopes1 - slopes2;
    pvalues(s_idx) = (1 + sum(slopes_diff < 0)) / (1 + nboot);
    
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
    bar(edges(1:end-1), hist1, .5, 'r', 'EdgeColor', 'none');
    bar(edges(1:end-1)+.5*edgediff, hist2, .5, 'b', 'EdgeColor', 'none');
    xlabel('PK slope');
    ylabel('bootstrap histogram');
    set(gca, 'YTick', 0);
    xlim([-.2, .2]);
    
    if s_idx == 1
        [~, label1] = getStairVar(phases(1));
        [~, label2] = getStairVar(phases(2));
        legend(label1, label2);
    end
end

end